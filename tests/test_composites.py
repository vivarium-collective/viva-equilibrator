"""Tests for the eQuilibrator composite generator."""

import pytest
from process_bigraph import Composite, allocate_core, gather_emitter_results

import viva_equilibrator  # noqa: F401  (ensure discovery registers the package)
from viva_equilibrator.composites.feasibility import equilibrator_feasibility_gate


def _engine_available() -> bool:
    try:
        from viva_equilibrator.processes import get_component_contribution
        get_component_contribution()
        return True
    except Exception:
        return False


def test_generator_is_registered():
    """The composite generator is discoverable in the shared registry."""
    from process_bigraph.composite_generator import _REGISTRY
    matches = [eid for eid in _REGISTRY if eid.endswith(".equilibrator_feasibility_gate")]
    assert matches, f"equilibrator_feasibility_gate missing; have {list(_REGISTRY)[:5]}"


def test_step_auto_registered():
    """EquilibratorThermoStep is auto-discovered (no manual register_link)."""
    core = allocate_core()
    assert "EquilibratorThermoStep" in core.link_registry


@pytest.mark.skipif(not _engine_available(), reason="equilibrator_api cache unavailable")
def test_feasibility_gate_composite_runs_and_emits():
    """The composite builds, runs, and emits exergonic ΔG' at vent conditions."""
    core = allocate_core()
    doc = equilibrator_feasibility_gate(core, pH=9.5, temperature=333.15, H2=5e-3, CO2=2e-3)
    sim = Composite({"state": doc}, core=core)
    sim.run(3.0)

    rows = gather_emitter_results(sim)[("emitter",)]
    assert rows, "no emitter rows"
    last = rows[-1]
    for name in ("formate", "acetate", "pyruvate"):
        assert last[f"dg_{name}"] < 0.0
        assert 0.0 <= last[f"feasible_{name}"] <= 1.0
