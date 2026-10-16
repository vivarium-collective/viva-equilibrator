"""Composite generators for the eQuilibrator wrapper."""

from . import feasibility  # noqa: F401  (side-effect: registers the generator)

from .feasibility import equilibrator_feasibility_gate

__all__ = ["equilibrator_feasibility_gate"]
