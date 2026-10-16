"""eQuilibrator process-bigraph wrapper.

A genuine bridge to the **real** ``equilibrator_api`` component-contribution
engine (with the offline ``equilibrator-cache-data`` database — no network).
Two public surfaces:

* :class:`EquilibratorThermoStep` — a per-tick ``Step`` that, given fluid
  conditions (pH, ionic strength, temperature) and the live concentrations of
  a configured reaction set's substrates/products, emits each reaction's
  transformed Gibbs energy ΔG' = ΔG'°(pH, I, T) + RT·ln(Q) and a feasibility
  multiplier in [0, 1] (logistic in −ΔG'/RT). ΔG'° is cached and recomputed
  only when pH/I/T change, so the Step is cheap to run every tick.
* :func:`standard_dg_primes` — a build-time helper that returns the real
  transformed ΔG'° (value + uncertainty) for a set of reaction formulas at
  given conditions, for reaction-universe / thermodynamic-table generation.

The ``equilibrator_api`` import is **lazy** (inside the engine accessor), so the
module imports and the Step registers even in an environment where the heavy
cache has not been materialised; only *running* the real engine touches it.
"""

from __future__ import annotations

import math

from process_bigraph import Step

# --------------------------------------------------------------------------- #
# Port schema fragments
# --------------------------------------------------------------------------- #
_MOLAR = {"_type": "float", "_units": "molar"}
_K = {"_type": "float", "_units": "kelvin"}
_DIMLESS = {"_type": "float", "_units": "dimensionless"}
# Transformed Gibbs energy sensor (overwrite): kilojoule/mole.
_KJ = {"_type": "overwrite[float]", "_units": "kilojoule/mole"}
_OW_DIMLESS = {"_type": "overwrite[float]", "_units": "dimensionless"}

_R = 8.314462618e-3  # gas constant, kJ/mol/K

# --------------------------------------------------------------------------- #
# Default reaction set: H2-driven reductive CO2 fixation at an alkaline vent.
# --------------------------------------------------------------------------- #
# Each entry: formula (eQuilibrator reaction string, KEGG ids) + the substrate
# stoichiometry used to form the reaction quotient Q = [product] / prod([sub]^n)
# + the product species whose concentration enters Q. KEGG ids: CO2 C00011,
# H2 C00282, H2O C00001, formate C00058, acetate C00033, pyruvate C00022.
VENT_CARBON_FIXATION = {
    "formate": {
        "formula": "kegg:C00011 + kegg:C00282 = kegg:C00058",
        "substrates": {"CO2": 1, "H2": 1},
        "product": "formate",
    },
    "acetate": {
        "formula": "2 kegg:C00011 + 4 kegg:C00282 = kegg:C00033 + 2 kegg:C00001",
        "substrates": {"CO2": 2, "H2": 4},
        "product": "acetate",
    },
    "pyruvate": {
        "formula": "3 kegg:C00011 + 5 kegg:C00282 = kegg:C00022 + 3 kegg:C00001",
        "substrates": {"CO2": 3, "H2": 5},
        "product": "pyruvate",
    },
}


# --------------------------------------------------------------------------- #
# Engine accessor + build-time helper
# --------------------------------------------------------------------------- #
_CC = None


def get_component_contribution():
    """Return a process-wide cached real ``ComponentContribution`` engine.

    Lazily constructs (and caches) the engine on first call; the heavy offline
    ``equilibrator-cache-data`` database is loaded once and reused. Raises
    ``ImportError`` if ``equilibrator_api`` is not installed.
    """
    global _CC
    if _CC is None:
        from equilibrator_api import ComponentContribution  # lazy, offline
        _CC = ComponentContribution()
    return _CC


def standard_dg_primes(formulas, p_h=7.0, ionic_strength_M=0.25,
                       temperature_K=298.15):
    """Real transformed ΔG'° for a set of reactions at ``(pH, I, T)``.

    Parameters
    ----------
    formulas : Mapping[str, str]
        ``{name: eQuilibrator reaction formula}`` (KEGG/compound ids, e.g.
        ``"kegg:C00011 + kegg:C00282 = kegg:C00058"``).
    p_h, ionic_strength_M, temperature_K : float
        Physiological conditions for the transformed (ΔG'°) estimate.

    Returns
    -------
    dict[str, tuple[float, float]]
        ``{name: (dg_prime_kj_mol, sigma_kj_mol)}`` — the component-contribution
        value and its uncertainty, both in kJ/mol.
    """
    from equilibrator_api import Q_
    cc = get_component_contribution()
    cc.p_h = Q_(float(p_h))
    cc.ionic_strength = Q_(f"{max(float(ionic_strength_M), 1e-6)} M")
    cc.temperature = Q_(f"{float(temperature_K)} K")
    out = {}
    for name, formula in formulas.items():
        rxn = cc.parse_reaction_formula(formula)
        dg = cc.standard_dg_prime(rxn)
        out[name] = (
            float(dg.value.to("kJ/mol").magnitude),
            float(dg.error.to("kJ/mol").magnitude),
        )
    return out


# --------------------------------------------------------------------------- #
# The per-tick feasibility-gate Step
# --------------------------------------------------------------------------- #
class EquilibratorThermoStep(Step):
    """Real eQuilibrator feasibility gate over a configured reaction set.

    For each reaction it emits ``dg_<name>`` (ΔG', kJ/mol) and
    ``feasible_<name>`` (a multiplier in [0, 1], ~1 when ΔG' is very negative
    and ~0 when positive), plus ``free_energy_logged`` as a witness that the
    thermodynamics were evaluated this tick. A ``Step`` (stateless per tick):
    it reads the conditions and reports sensor values. ΔG'° is cached per
    rounded ``(pH, I, T)`` so only ΔG''s ln(Q) term is recomputed each tick.

    The reaction set defaults to :data:`VENT_CARBON_FIXATION` (H2-driven
    reductive CO2 fixation to formate/acetate/pyruvate); pass ``reactions`` in
    config to gate a different network. Input/output ports are derived from the
    reaction set, so the default config reproduces the ports
    ``{pH, ionic_strength, temperature, CO2, H2, formate, acetate, pyruvate}``
    → ``{dg_*, feasible_*, free_energy_logged}``.
    """

    contract = {
        "summary": (
            "Real eQuilibrator feasibility gate: (pH, I, T, substrates, products) "
            "-> dG' + feasibility per reaction."
        ),
        "description": (
            "A stateless Step bridging the real equilibrator_api component-contribution "
            "engine (offline cache). Each tick it reads the fluid conditions and the live "
            "concentrations of the configured reaction set and emits the transformed Gibbs "
            "energy dG' = dG'0(pH, I, T) + RT ln Q plus a logistic feasibility multiplier "
            "in [0,1] that downstream kinetics can use as a thermodynamic gate. dG'0 is "
            "cached and recomputed only when pH / ionic strength / temperature change."
        ),
        "math": ["dG' = dG'0 + RT ln Q", "feasible = 1 / (1 + exp(dG'/RT))"],
        "symbols": {"Q": "reaction quotient", "R": "gas constant", "T": "temperature (K)"},
        "assumptions": [
            "Physiological (transformed) dG' at the given pH / I / T.",
            "Q = [product] / prod([substrate]^stoich); water activity ~ 1.",
        ],
        "references": ["noor2013consistent", "beber2022equilibrator", "preiner2020hydrogen"],
    }

    config_schema = {
        # Activity floor for ln(Q) so a zero concentration can't blow up the log.
        "conc_floor": {"_type": "float", "_default": 1.0e-9, "_units": "molar"},
        # Reaction set: {name: {formula, substrates:{species:stoich}, product}}.
        # ``quote`` keeps the nested dict out of the engine's type machinery.
        # The default is empty and read as "use VENT_CARBON_FIXATION" — a
        # non-empty default would be *merged* into (not replaced by) a provided
        # reaction set, silently re-adding the vent reactions.
        "reactions": {"_type": "quote", "_default": {}},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        reactions = self.config.get("reactions") or VENT_CARBON_FIXATION
        self._reactions = dict(reactions)
        # Union of substrate species and the set of product species.
        subs = set()
        for spec in self._reactions.values():
            subs.update(spec["substrates"].keys())
        self._substrates = sorted(subs)
        self._products = [spec["product"] for spec in self._reactions.values()]
        self._parsed = None               # {name: eQuilibrator reaction}
        self._dg0_cache = {}              # {(pH,I,T rounded): {name: dg0 kJ/mol}}

    def inputs(self):
        ports = {"pH": _DIMLESS, "ionic_strength": _MOLAR, "temperature": _K}
        for species in self._substrates:
            ports[species] = _MOLAR
        for species in self._products:
            ports[species] = _MOLAR
        return ports

    def outputs(self):
        ports = {}
        for name in self._reactions:
            ports[f"dg_{name}"] = _KJ
            ports[f"feasible_{name}"] = _OW_DIMLESS
        ports["free_energy_logged"] = {"_type": "overwrite[boolean]"}
        return ports

    def _parse(self):
        if self._parsed is None:
            cc = get_component_contribution()
            self._parsed = {
                name: cc.parse_reaction_formula(spec["formula"])
                for name, spec in self._reactions.items()
            }
        return self._parsed

    def _dg0(self, ph: float, ionic: float, temp_k: float) -> dict:
        key = (round(ph, 2), round(ionic, 3), round(temp_k, 1))
        cached = self._dg0_cache.get(key)
        if cached is not None:
            return cached
        from equilibrator_api import Q_
        cc = get_component_contribution()
        parsed = self._parse()
        cc.p_h = Q_(float(key[0]))
        cc.ionic_strength = Q_(f"{max(key[1], 1e-6)} M")
        cc.temperature = Q_(f"{key[2]} K")
        out = {
            name: float(cc.standard_dg_prime(rxn).value.to("kJ/mol").magnitude)
            for name, rxn in parsed.items()
        }
        self._dg0_cache[key] = out
        return out

    def feasibility(self, dg0: dict, conc: dict, temp_k: float) -> dict:
        """Pure thermodynamic math: ``{name: (dg, feasible)}`` from cached ΔG'°.

        Separated from the engine so it is unit-testable without the cache.
        ``conc`` maps every substrate and product species to a floored molar
        concentration.
        """
        rt = _R * temp_k
        out = {}
        for name, spec in self._reactions.items():
            denom = 1.0
            for species, stoich in spec["substrates"].items():
                denom *= conc[species] ** stoich
            q = conc[spec["product"]] / denom
            dg = dg0[name] + rt * math.log(q)
            feasible = 1.0 / (1.0 + math.exp(dg / rt))
            out[name] = (dg, feasible)
        return out

    def update(self, state):
        ph = float(state.get("pH", 7.0))
        ionic = float(state.get("ionic_strength", 0.25))
        temp_k = max(float(state.get("temperature", 298.15)), 250.0)
        floor = self.config["conc_floor"]
        conc = {
            species: max(float(state.get(species, 0.0)), floor)
            for species in set(self._substrates) | set(self._products)
        }

        dg0 = self._dg0(ph, ionic, temp_k)
        results = self.feasibility(dg0, conc, temp_k)

        out = {"free_energy_logged": True}
        for name, (dg, feasible) in results.items():
            out[f"dg_{name}"] = dg
            out[f"feasible_{name}"] = feasible
        return out
