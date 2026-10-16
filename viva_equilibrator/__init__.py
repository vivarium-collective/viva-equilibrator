"""viva-equilibrator: process-bigraph wrapper for eQuilibrator.

Bridges the real ``equilibrator_api`` component-contribution engine (offline
``equilibrator-cache-data``) into process-bigraph: a per-tick transformed-Gibbs
feasibility gate (:class:`EquilibratorThermoStep`) and a build-time ΔG'°
annotation helper (:func:`standard_dg_primes`).
"""

from viva_equilibrator.processes import (
    EquilibratorThermoStep,
    VENT_CARBON_FIXATION,
    standard_dg_primes,
    get_component_contribution,
)
from viva_equilibrator.types import register_equilibrator_types

__all__ = [
    "EquilibratorThermoStep",
    "VENT_CARBON_FIXATION",
    "standard_dg_primes",
    "get_component_contribution",
    "register_equilibrator_types",
]
