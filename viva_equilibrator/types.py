"""Custom bigraph-schema types for the eQuilibrator wrapper.

The wrapper uses built-in types only (unit-bearing ``float``/``overwrite``,
``quote`` for the reaction-set config). This hook exists for future custom
type registrations — e.g. a reaction-spec type or a ΔG'-with-uncertainty type.
"""


def register_equilibrator_types(core):
    """Register custom types used by the eQuilibrator processes. No-op for now."""
    return core
