"""Module containing the ``Species`` class and sub-classes.

Species
-------
In fluid simulations, the species are the chemical compounds present.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from ansys.dpf.post.simulation import Simulation


class Species:
    """In fluid simulation, defines a chemical species."""

    def __init__(self, name: str, id: int):
        """Initialize this class."""
        self._name = name
        self._id = id

    @property
    def name(self) -> str:
        """Return the name of the species."""
        return self._name

    @property
    def id(self) -> int:
        """Return the ID of the species."""
        return self._id

    def __repr__(self) -> str:
        """String representation of the instance."""
        return f'Species<name: "{self._name}", id={self._id}>'


class SpeciesList:
    """In fluid simulation, defines a list of chemical species."""

    def __init__(self, simulation: Simulation):
        """Initialize this class."""
        self._species = []
        if "species" in simulation.result_info.available_qualifier_labels:
            support = simulation.result_info.qualifier_label_support("species")
            species_names_field = support.string_field_support_by_property("names")
            for id, name in enumerate(species_names_field.data_as_list):
                self._species.append(Species(name, id))

    def __repr__(self) -> str:
        """String representation of the instance."""
        text = "["
        for species in self:
            text += repr(species) + ", "
        text += "]"
        return text

    def __len__(self):
        """Length of the instance."""
        return len(self._species)

    def __str__(self) -> str:
        """String representation of the instance."""
        text = f"{len(self)} species available\n"
        for species in self._species:
            text += f"{species.id}: {species.name}\n"
        return text

    def __getitem__(self, item: int) -> Species:
        """Return the Species at the given position in the list."""
        return self._species[item]
