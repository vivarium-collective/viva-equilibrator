# viva-equilibrator

Process-bigraph wrapper for **[eQuilibrator](https://equilibrator.weizmann.ac.il/)** —
the component-contribution thermodynamics engine (Noor 2013, Beber 2022). It
bridges the real `equilibrator_api` (with the offline `equilibrator-cache-data`
database — no network) into process-bigraph so transformed Gibbs energies and
reaction feasibility can drive a simulation.

## What it provides

| Surface | Kind | Use |
|---|---|---|
| `EquilibratorThermoStep` | `Step` | Per-tick ΔG′ feasibility gate over a configurable reaction set |
| `standard_dg_primes(...)` | function | Build-time ΔG′° (value + uncertainty) for reaction-table / universe generation |
| `get_component_contribution()` | function | Cached real `ComponentContribution` engine |
| `equilibrator_feasibility_gate` | composite generator | Demo: a fixed vent condition driving the gate |

### `EquilibratorThermoStep`

Given fluid conditions (pH, ionic strength, temperature) and the live
concentrations of a reaction set's substrates/products, it emits per reaction:

- `dg_<name>` — the transformed Gibbs energy ΔG′ = ΔG′°(pH, I, T) + RT·ln(Q), kJ/mol
- `feasible_<name>` — a multiplier in [0, 1], logistic in −ΔG′/RT (≈1 strongly exergonic, ≈0 endergonic)

plus `free_energy_logged` as a per-tick witness. ΔG′° is cached per rounded
`(pH, I, T)`, so only the ln(Q) term is recomputed each tick.

The reaction set defaults to **`VENT_CARBON_FIXATION`** — H₂-driven reductive CO₂
fixation to formate/acetate/pyruvate (the alkaline-hydrothermal-vent set). Pass a
`reactions` config (`{name: {formula, substrates: {species: stoich}, product}}`)
to gate any other network; the input/output ports derive from that set.

## Installation

```bash
# editable, for development
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

`equilibrator-api` downloads its offline cache once on first use; thereafter the
engine runs with no network. Once installed, the process registers automatically
via `bigraph_schema.package.discover` — no manual `register_link()` needed.

## Quick start

```python
from process_bigraph import Composite, allocate_core, gather_emitter_results
from viva_equilibrator.composites.feasibility import equilibrator_feasibility_gate

core = allocate_core()
doc = equilibrator_feasibility_gate(core, pH=9.5, temperature=333.15, H2=5e-3, CO2=2e-3)
sim = Composite({"state": doc}, core=core)
sim.run(3.0)
last = gather_emitter_results(sim)[("emitter",)][-1]
print(last["dg_formate"], last["feasible_formate"])   # e.g. -48.7 kJ/mol, 1.0
```

Build-time ΔG′° table:

```python
from viva_equilibrator import standard_dg_primes
dg = standard_dg_primes(
    {"formate": "kegg:C00011 + kegg:C00282 = kegg:C00058"},
    p_h=9.5, ionic_strength_M=0.25, temperature_K=333.15,
)
# {"formate": (value_kJ_per_mol, sigma_kJ_per_mol)}
```

## Used by

[`viva-genesis`](https://github.com/vivarium-collective/viva-genesis) imports this
package for its origins-of-life vent campaign: `EquilibratorThermoStep` is the
thermodynamic feasibility gate on carbon fixation, and `standard_dg_primes`
annotates the generated reaction universe with transformed Gibbs energies.

## References

- Noor et al. (2013), *Consistent estimation of Gibbs energy using component contributions*. PLoS Comput Biol.
- Beber et al. (2022), *eQuilibrator 3.0: a database solution for thermodynamic constant estimation*. Nucleic Acids Res.

## License

MIT
