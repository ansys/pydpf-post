# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Module containing the ``Species`` class and sub-classes.

Species
-------
In fluid simulations, the species are the chemical compounds present.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:  # pragma: no cover
    from ansys.dpf.post.simulation import Simulation


class Species:
    """In fluid simulation, defines a chemical species."""

    def __init__(self, name: str, id: int):  # pylint: disable=redefined-builtin
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


class SpeciesDict:
    """Dictionary of species available in the fluid simulation.

    Accepts either species name or species ID as key.
    """

    def __init__(self, simulation: Simulation):
        """Initialize this class."""
        self._species = []
        self._ids = []
        self._idx = 0
        self._names = []
        if "species" in simulation.result_info.available_qualifier_labels:
            support = simulation.result_info.qualifier_label_support("species")
            species_names_field = support.string_field_support_by_property("names")
            names = species_names_field.data_as_list
            ids = species_names_field.scoping.ids
            for i, name in enumerate(names):
                self._species.append(Species(name, ids[i]))
                self._names.append(name)
                self._ids.append(ids[i])

    def __repr__(self) -> str:
        """String representation of the instance."""
        text = "{"
        for species in self:
            text += repr(species) + ", "
        text += "}"
        return text

    def __len__(self):
        """Length of the instance."""
        return len(self._species)

    def __str__(self) -> str:
        """String representation of the instance."""
        text = f"{len(self)} species available\n"
        text += repr(self)
        return text

    def __getitem__(self, item: Union[int, str]) -> Species:
        """Returns the Species with the given name or ID."""
        try:
            if isinstance(item, str):
                index = self._names.index(item)
            else:
                index = self._ids.index(item)
            return self._species[index]
        except Exception:
            raise ValueError(f"{item} is not a valid Species ID or Species name.")

    def __next__(self) -> Species:
        """Returns the next entity in the iterator."""
        if self._idx >= len(self):
            raise StopIteration
        out = self[self._ids[self._idx]]
        self._idx += 1
        return out

    def __iter__(self) -> SpeciesDict:
        """Returns the object to iterate over."""
        self._idx = 0
        return self
