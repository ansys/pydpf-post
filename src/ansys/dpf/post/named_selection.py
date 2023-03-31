"""This module contains NamedSelectionsDict, NamedSelectionsIterator and NamedSelection classes."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import List

import ansys.dpf.core as dpf


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

    # def __setitem__(self, key: str, value: List[int]):
    #    """Implements [] setter access function."""
    #    self._meshed_region.set_named_selection_scoping(
    #        named_selection_name=key, scoping=dpf.Scoping(ids=value)
    #    )

    def __len__(self) -> int:
        """Returns the length of the dictionary (number of named selections)."""
        return len(self.keys())

    def keys(self) -> List[str]:
        """Returns the available named selections."""
        return self._meshed_region.available_named_selections

    def values(self) -> List[NamedSelection]:
        """Returns list of the values of all the named selections."""
        return [self[key] for key in self.keys()]

    def has_key(self, key) -> bool:
        """Returns True the given key is present in available named selections."""
        return key in self.keys()

    def __delitem__(self, __key):
        """Not implemented."""
        pass

    def __iter__(self) -> NamedSelectionsIterator:
        """Returns an iterator to access this dictionary."""
        return NamedSelectionsIterator(self)


class NamedSelection(dpf.Scoping):
    """Class decorating dpf.Scoping with a name attribute."""

    def __init__(self, name: str, scoping: dpf.Scoping):
        """Constructs a NamedSelection from a name and a Scoping."""
        super().__init__(scoping)

        self._name = name

    @property
    def name(self) -> str:
        """Returns the name."""
        return self._name

    # @name.setter
    # def name(self, val: str):
    #    self._name = val

    def __repr__(self) -> str:
        """Pretty print string of the NamedSelection."""
        return f"NamedSelection '{self.name}' with scoping {repr(super())}"
