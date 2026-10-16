"""Tests for the eQuilibrator process-bigraph wrapper.

The pure-math tests run everywhere. The real-engine tests are guarded: they
skip (not fail) when ``equilibrator_api`` or its offline cache is unavailable,
so CI without the heavy cache stays green while a provisioned machine exercises
the genuine component-contribution path.
"""

import math

import pytest
from process_bigraph import allocate_core

from viva_equilibrator.processes import (
    EquilibratorThermoStep,
    VENT_CARBON_FIXATION,
    standard_dg_primes,
)


def _engine_available() -> bool:
    try:
        from viva_equilibrator.processes import get_component_contribution
        get_component_contribution()
        return True
    except Exception:
        return False


def test_default_ports_match_vent_carbon_fixation():
    """Default reaction set yields the vent carbon-fixation port surface."""
    core = allocate_core()
    step = EquilibratorThermoStep(config={}, core=core)

    ins = set(step.inputs())
    outs = set(step.outputs())

    assert ins == {"pH", "ionic_strength", "temperature",
                   "CO2", "H2", "formate", "acetate", "pyruvate"}
    assert outs == {"dg_formate", "dg_acetate", "dg_pyruvate",
                    "feasible_formate", "feasible_acetate", "feasible_pyruvate",
                    "free_energy_logged"}


def test_custom_reaction_set_drives_ports():
    """A one-reaction config produces exactly that reaction's ports."""
    core = allocate_core()
    rxns = {"formate": VENT_CARBON_FIXATION["formate"]}
    step = EquilibratorThermoStep(config={"reactions": rxns}, core=core)

    assert set(step.inputs()) == {"pH", "ionic_strength", "temperature", "CO2", "H2", "formate"}
    assert set(step.outputs()) == {"dg_formate", "feasible_formate", "free_energy_logged"}


def test_feasibility_logistic_is_pure_math():
    """feasibility() is a logistic in -dG'/RT with no engine dependency."""
    core = allocate_core()
    step = EquilibratorThermoStep(config={}, core=core)

    temp_k = 333.15
    conc = {"CO2": 2e-3, "H2": 5e-3, "formate": 1e-6, "acetate": 1e-6, "pyruvate": 1e-6}
    # Very negative dG'0 -> feasible ~ 1; very positive -> feasible ~ 0.
    dg0 = {"formate": -50.0, "acetate": +50.0, "pyruvate": 0.0}
    results = step.feasibility(dg0, conc, temp_k)

    assert results["formate"][1] > 0.99
    assert results["acetate"][1] < 0.01
    # At dG' == 0 (pyruvate dg0=0 and ln Q term), feasibility is the logistic midpoint-ish;
    # just assert it is a valid probability.
    for _dg, feasible in results.values():
        assert 0.0 <= feasible <= 1.0


def test_feasibility_matches_closed_form():
    """feasibility() equals dg0 + RT ln Q and the logistic of it."""
    core = allocate_core()
    step = EquilibratorThermoStep(config={"reactions": {"formate": VENT_CARBON_FIXATION["formate"]}},
                                  core=core)
    temp_k = 300.0
    R = 8.314462618e-3
    rt = R * temp_k
    conc = {"CO2": 1e-2, "H2": 1e-2, "formate": 1e-4}
    dg0 = {"formate": -30.0}
    q = conc["formate"] / (conc["CO2"] ** 1 * conc["H2"] ** 1)
    expected_dg = dg0["formate"] + rt * math.log(q)
    expected_feasible = 1.0 / (1.0 + math.exp(expected_dg / rt))

    dg, feasible = step.feasibility(dg0, conc, temp_k)["formate"]
    assert abs(dg - expected_dg) < 1e-9
    assert abs(feasible - expected_feasible) < 1e-9


@pytest.mark.skipif(not _engine_available(), reason="equilibrator_api cache unavailable")
def test_real_engine_update_is_exergonic_in_h2_rich_vent():
    """Real eQuilibrator: H2-rich reductive CO2 fixation is strongly exergonic."""
    core = allocate_core()
    step = EquilibratorThermoStep(config={}, core=core)
    out = step.update({
        "pH": 9.5, "ionic_strength": 0.25, "temperature": 333.15,
        "CO2": 2e-3, "H2": 5e-3, "formate": 1e-6, "acetate": 1e-6, "pyruvate": 1e-6,
    })
    assert out["free_energy_logged"] is True
    for name in ("formate", "acetate", "pyruvate"):
        assert out[f"dg_{name}"] < 0.0                 # exergonic
        assert 0.0 <= out[f"feasible_{name}"] <= 1.0


@pytest.mark.skipif(not _engine_available(), reason="equilibrator_api cache unavailable")
def test_standard_dg_primes_returns_value_and_uncertainty():
    """standard_dg_primes returns (value, sigma) kJ/mol for each formula."""
    out = standard_dg_primes(
        {"formate": VENT_CARBON_FIXATION["formate"]["formula"]},
        p_h=7.0, ionic_strength_M=0.25, temperature_K=298.15,
    )
    value, sigma = out["formate"]
    assert isinstance(value, float) and isinstance(sigma, float)
    assert sigma >= 0.0
    assert -100.0 < value < 100.0                      # sane magnitude


@pytest.mark.skipif(not _engine_available(), reason="equilibrator_api cache unavailable")
def test_dg0_is_cached_per_condition():
    """dG'0 is recomputed only when (pH, I, T) changes."""
    core = allocate_core()
    step = EquilibratorThermoStep(config={}, core=core)
    step.update({"pH": 9.5, "ionic_strength": 0.25, "temperature": 333.15,
                 "CO2": 2e-3, "H2": 5e-3})
    assert len(step._dg0_cache) == 1
    step.update({"pH": 9.5, "ionic_strength": 0.25, "temperature": 333.15,
                 "CO2": 1e-3, "H2": 9e-3})             # same conditions, new conc
    assert len(step._dg0_cache) == 1                   # no recompute
    step.update({"pH": 7.0, "ionic_strength": 0.25, "temperature": 333.15,
                 "CO2": 2e-3, "H2": 5e-3})             # new pH
    assert len(step._dg0_cache) == 2
