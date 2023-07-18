"""This module contains NamedSelectionsDict, NamedSelectionsIterator and NamedSelection classes."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, List, Union

import ansys.dpf.core as dpf
import ansys.dpf.core.dpf_array as dpf_array

if TYPE_CHECKING:  # pragma: no cover
    from ansys.dpf.post.mesh import Mesh


class NamedSelections(Mapping):
    """Dictionary of available named selections for a given mesh."""

    def __init__(self, mesh: Mesh):
        """Initialize Named Selections dictionary from a Mesh."""
        self._idx = 0
        self._meshed_region = mesh._meshed_region

    def __getitem__(self, key: str) -> NamedSelection:
        """Get named selection of name equal to the given key."""
        if key in self._meshed_region.available_named_selections:
            scoping = self._meshed_region.named_selection(key)
            return NamedSelection(key, scoping)

        raise KeyError(f'named selection "{key}" could not be found')

    def keys(self) -> List[str]:
        """Returns the available named selections."""
        return self._meshed_region.available_named_selections

    def __len__(self) -> int:
        """Returns the length of the dictionary (number of named selections)."""
        return len(self.keys())

    def __iter__(self) -> NamedSelections:
        """Get base iterator."""
        self._idx = 0
        return self

    def __next__(self) -> str:
        """Returns next value."""
        if self._idx >= len(self):
            raise StopIteration
        res = self.keys()[self._idx]
        self._idx += 1
        return res

    def __str__(self) -> str:
        """String representation."""
        txt = f"NamedSelections dictionary with {len(self)} named selections:"
        for ns in self:
            txt += f"\n\t- '{ns}'"
        return txt


class NamedSelection:
    """Named Selection class associating a name to a list of mesh entities."""

    def __init__(
        self,
        name: str,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Constructs a NamedSelection from a name and a list of entity IDs."""
        tot = (
            (node_ids is not None)
            + (element_ids is not None)
            + (face_ids is not None)
            + (cell_ids is not None)
            + (scoping is not None)
        )
        if tot > 1:
            raise ValueError(
                "NamedSelection accepts only one argument giving a list of entities."
            )
        elif tot == 0:
            raise ValueError("No list of entity IDs given.")
        if scoping:
            self._scoping = scoping
        elif node_ids:
            self._scoping = dpf.mesh_scoping_factory.nodal_scoping(node_ids=node_ids)
        elif element_ids:
            self._scoping = dpf.mesh_scoping_factory.elemental_scoping(
                element_ids=element_ids
            )
        elif face_ids:
            self._scoping = dpf.mesh_scoping_factory.face_scoping(face_ids=face_ids)
        elif cell_ids:
            self._scoping = dpf.mesh_scoping_factory.elemental_scoping(
                element_ids=cell_ids
            )

        self._name = name

    def __eq__(self, other) -> bool:
        """Tests equality of location and IDs."""
        return (self.location == other.location) and (all(self.ids == other.ids))

    @property
    def name(self) -> str:
        """Returns the name."""
        return self._name

    # Scoping forwarding
    def set_id(self, index: int, scopingid: int):
        """Sets the ID of the underlying scoping's index."""
        self._scoping.set_id(index, scopingid)

    def id(self, index: int) -> int:
        """Retrieve the ID at a given index."""
        return self._scoping.id(index)

    def index(self, id: int) -> int:  # pylint: disable=redefined-builtin
        """Retrieve the index of a given ID."""
        return self._scoping.index(id)

    @property
    def ids(self) -> Union[dpf_array.DPFArray, List[int]]:
        """Retrieve a list of IDs in the underlying scoping."""
        return self._scoping.ids

    @property
    def location(self) -> str:
        """Location of the IDs as a string."""
        return self._scoping.location

    @property
    def size(self) -> int:
        """Length of the list of IDs."""
        return self._scoping.size

    def __repr__(self) -> str:
        """Pretty print string of the NamedSelection."""
        return f"NamedSelection '{self.name}'\n with {self._scoping.__str__()}"
