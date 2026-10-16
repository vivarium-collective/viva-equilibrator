"""A minimal composite: a fixed vent environment feeding the ΔG' feasibility gate.

Holds the fluid conditions and dissolved-pool concentrations in a ``world``
store, runs :class:`~viva_equilibrator.processes.EquilibratorThermoStep` over
them, and emits the per-reaction ΔG' and feasibility each tick. A read-only
demonstration of the real eQuilibrator engine inside process-bigraph.
"""

from process_bigraph.composite_generator import composite_generator


_EMIT = [
    "gibbs/dg_formate", "gibbs/dg_acetate", "gibbs/dg_pyruvate",
    "gibbs/feasible_formate", "gibbs/feasible_acetate", "gibbs/feasible_pyruvate",
    "global_time",
]


@composite_generator(
    name="equilibrator_feasibility_gate",
    description=(
        "Real eQuilibrator ΔG' feasibility gate over H2-driven reductive CO2 "
        "fixation (formate/acetate/pyruvate) at a fixed vent condition."
    ),
    parameters={
        "pH": {"type": "float", "default": 9.5, "description": "Fluid pH"},
        "temperature": {"type": "float", "default": 333.15,
                        "description": "Temperature (K)"},
        "ionic_strength": {"type": "float", "default": 0.25,
                           "description": "Ionic strength (molar)"},
        "H2": {"type": "float", "default": 5.0e-3,
               "description": "Dissolved H2 (molar)"},
        "CO2": {"type": "float", "default": 2.0e-3,
                "description": "Dissolved CO2 (molar)"},
    },
    emitters=[{"address": "local:ParquetEmitter", "paths": _EMIT}],
)
def equilibrator_feasibility_gate(core=None, *, pH=9.5, temperature=333.15,
                                  ionic_strength=0.25, H2=5.0e-3, CO2=2.0e-3):
    return {
        "thermo": {
            "_type": "step",
            "address": "local:EquilibratorThermoStep",
            "config": {},
            "inputs": {
                "pH": ["world", "pH"],
                "ionic_strength": ["world", "ionic_strength"],
                "temperature": ["world", "temperature"],
                "CO2": ["world", "CO2"],
                "H2": ["world", "H2"],
                "formate": ["world", "formate"],
                "acetate": ["world", "acetate"],
                "pyruvate": ["world", "pyruvate"],
            },
            "outputs": {
                "dg_formate": ["gibbs", "dg_formate"],
                "dg_acetate": ["gibbs", "dg_acetate"],
                "dg_pyruvate": ["gibbs", "dg_pyruvate"],
                "feasible_formate": ["gibbs", "feasible_formate"],
                "feasible_acetate": ["gibbs", "feasible_acetate"],
                "feasible_pyruvate": ["gibbs", "feasible_pyruvate"],
                "free_energy_logged": ["gibbs", "free_energy_logged"],
            },
        },
        "world": {
            "pH": pH, "ionic_strength": ionic_strength, "temperature": temperature,
            "H2": H2, "CO2": CO2,
            "formate": 1.0e-6, "acetate": 1.0e-6, "pyruvate": 1.0e-6,
        },
        "gibbs": {
            "dg_formate": 0.0, "dg_acetate": 0.0, "dg_pyruvate": 0.0,
            "feasible_formate": 0.0, "feasible_acetate": 0.0, "feasible_pyruvate": 0.0,
            "free_energy_logged": False,
        },
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {
                "dg_formate": "float", "dg_acetate": "float", "dg_pyruvate": "float",
                "feasible_formate": "float", "feasible_acetate": "float",
                "feasible_pyruvate": "float", "global_time": "float",
            }},
            "inputs": {
                "dg_formate": ["gibbs", "dg_formate"],
                "dg_acetate": ["gibbs", "dg_acetate"],
                "dg_pyruvate": ["gibbs", "dg_pyruvate"],
                "feasible_formate": ["gibbs", "feasible_formate"],
                "feasible_acetate": ["gibbs", "feasible_acetate"],
                "feasible_pyruvate": ["gibbs", "feasible_pyruvate"],
                "global_time": ["global_time"],
            },
        },
    }
