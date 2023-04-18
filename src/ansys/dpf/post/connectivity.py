"""Module containing wrapper class for the connectivities property fields."""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum
from typing import List

from ansys.dpf.core.property_field import PropertyField
from ansys.dpf.core.scoping import Scoping


class Mode(Enum):
    """Enum made for internal use, to dictate the behavior of ConnectivityList."""

    IDS_FROM_IDX = 1
    IDS_FROM_ID = 2
    IDX_FROM_IDX = 3
    IDX_FROM_ID = 4


class ConnectivityList(Sequence):
    """Very basic wrapper around elemental and nodal connectivities fields."""

    def __init__(self, field: PropertyField, scoping: Scoping, mode: Mode):
        """Constructs a ConnectivityList by wrapping given PropertyField."""
        self._field = field
        self._mode = mode
        self._scoping = scoping

        self.local_scoping = None

    def __getitem__(self, key: int) -> List[int]:
        """Returns a list of indexes or IDs for a given index or ID, see Mode Enum."""
        if self._mode == Mode.IDS_FROM_IDX:
            return self._get_ids_from_idx(key)
        elif self._mode == Mode.IDS_FROM_ID:
            return self._get_ids_from_id(key)
        elif self._mode == Mode.IDX_FROM_IDX:
            return self._get_idx_from_idx(key)
        elif self._mode == Mode.IDX_FROM_ID:
            return self._get_idx_from_id(key)

    def _get_ids_from_idx(self, idx: int) -> List[int]:
        """Helper method to retrieve list of IDs from a given index."""
        return self._to_ids(self._field.get_entity_data(idx))

    def _get_ids_from_id(self, id: int) -> List[int]:
        """Helper method to retrieve list of IDs from a given ID."""
        return self._to_ids(self._field.get_entity_data_by_id(id))

    def _get_idx_from_idx(self, idx: int) -> List[int]:
        """Helper method to retrieve list of indexes from a given index."""
        return self._field.get_entity_data(idx)

    def _get_idx_from_id(self, id: int) -> List[int]:
        """Helper method to retrieve list of indexes from a given ID."""
        return self._field.get_entity_data_by_id(id)

    @property
    def by_id(self) -> ConnectivityList:
        """Returns an equivalent list which accepts an ID instead of an index in __getitem__."""
        if self._mode == Mode.IDS_FROM_IDX:
            return ConnectivityList(self._field, self._scoping, Mode.IDS_FROM_ID)
        elif self._mode == Mode.IDX_FROM_IDX:
            return ConnectivityList(self._field, self._scoping, Mode.IDX_FROM_ID)
        return self

    def _to_ids(self, indices) -> List[int]:
        """Helper method to convert a list of indexes into a list of IDs."""
        if not self.local_scoping:
            self.local_scoping = self._scoping.as_local_scoping()

        to_id = self.local_scoping.id
        return list(map(to_id, indices))

    def __len__(self) -> int:
        """Returns the number of entities."""
        return len(self._field._get_data_pointer())
