"""Module containing the ``Phase`` class and sub-classes.

Phase
-----
In fluid simulations, a phase defines the physical state of a given species.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ansys.dpf.post.simulation import Simulation


class Phase:
    """Physical state of a given species in a fluid simulation."""

    def __init__(self, name: str, id: int):
        """Initialize this class."""
        self._name = name
        self._id = id

    @property
    def name(self) -> str:
        """Return the name of the phase."""
        return self._name

    @property
    def id(self) -> int:
        """Return the ID of the phase."""
        return self._id

    def __repr__(self) -> str:
        """String representation of the instance."""
        return f'Phase<name: "{self._name}", id={self._id}>'


class Phases:
    """List of physical states in a fluid simulation."""

    def __init__(self, simulation: Simulation):
        """Initialize this class."""
        self._phases = []
        if "phase" in simulation.result_info.available_qualifier_labels:
            phase_support = simulation.result_info.qualifier_label_support("phase")
            phase_names_field = phase_support.string_field_support_by_property("names")
            for id, name in enumerate(phase_names_field.data_as_list):
                self._phases.append(Phase(name, id))

    def __repr__(self) -> str:
        """String representation of the instance."""
        text = "["
        for phase in self:
            text += repr(phase) + ", "
        text += "]"
        return text

    def __len__(self):
        """Length of the instance."""
        return len(self._phases)

    def __str__(self) -> str:
        """String representation of the instance."""
        text = f"{len(self)} phases available\n"
        for phase in self._phases:
            text += f"{phase.id}: {phase.name}\n"
        return text

    def __getitem__(self, item: int) -> Phase:
        """Return the Phase at the given position in the list."""
        return self._phases[item]
