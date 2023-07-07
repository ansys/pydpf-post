"""This module contains NamedSelectionsDict, NamedSelectionsIterator and NamedSelection classes."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
import copy
from typing import List, Union

import ansys.dpf.core as dpf
import ansys.dpf.core.dpf_array as dpf_array


class NamedSelectionsIterator(Iterator):
    """Iterator implementation for NamedSelectionsDict."""

    def __init__(self, ns_dict: NamedSelectionsDict):
        """Initialize the Named Selection Iterator. see NamedSelectionsDict."""
        self.idx = 0
        self.ns_dict = ns_dict

    def __iter__(self) -> NamedSelectionsIterator:
        """Get base iterator."""
        self.idx = 0
        return self

    def __next__(self) -> str:
        """Returns next value."""
        if self.idx < len(self.ns_dict):
            res = self.ns_dict.keys()[self.idx]
            self.idx += 1
            return res
        else:
            raise StopIteration


class NamedSelectionsDict(Mapping):
    """Proxy class to expose Named Selections interface to post.Mesh."""

    def __init__(self, meshed_region: dpf.MeshedRegion):
        """Initialize Named Selections dictionary from internal Meshed Region."""
        self._meshed_region = meshed_region

    def __getitem__(self, key: str) -> NamedSelection:
        """Implements [] getter access function."""
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

    def __delitem__(self, __key):
        """Not implemented."""
        return NotImplementedError

    def __iter__(self) -> NamedSelectionsIterator:
        """Returns an iterator to access this dictionary."""
        return NamedSelectionsIterator(self)


class NamedSelection:
    """Class decorating dpf.Scoping with a name attribute."""

    def __init__(self, name: str, scoping: dpf.Scoping):
        """Constructs a NamedSelection from a name and a Scoping."""
        self._scoping = scoping
        self._name = name

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

    def deep_copy(self, server=None) -> NamedSelection:
        """Create a deep copy of the underlying scoping's data on a given server."""
        new_scoping = self._scoping.deep_copy(server)
        new_name = copy.copy(self._name)
        return NamedSelection(new_name, new_scoping)

    def as_local_scoping(self) -> NamedSelection:
        """Create a deep copy of the underlying scoping that can be modified locally."""
        local_scoping = self._scoping.as_local_scoping()
        local_name = copy.copy(self._name)
        return NamedSelection(local_name, local_scoping)

    def __repr__(self) -> str:
        """Pretty print string of the NamedSelection."""
        return f"NamedSelection '{self.name}'\n with {self._scoping.__str__()}"
