"""Module containing wrapper class for the connectivities property fields."""

from __future__ import annotations

from collections.abc import Collection, Iterator
from enum import Enum
from typing import List

from ansys.dpf.core.property_field import PropertyField
from ansys.dpf.core.scoping import Scoping


class ReturnMode(Enum):
    """Enum made for internal use, to dictate the behavior of ConnectivityList."""
    IDS = 1
    IDX = 2

class ConnectivityListIterator(Iterator):
    def __init__(self, conn_list: ConnectivityList):
        self._conn_list = conn_list
        self._idx = 0
    
    def __next__(self) -> Any:
        if self._idx >= self._conn_list.__len__():
            raise StopIteration
        
        ret = self._conn_list[self._idx]
        self._idx += 1
        return ret
    
    def __iter__(self) -> ConnectivityListIterator:
        return ConnectivityListIdx(self._conn_list)

class ConnectivityListIdx(Collection):
    """Very basic wrapper around elemental and nodal connectivities fields."""

    def __init__(self, field: PropertyField, scoping: Scoping, mode: ReturnMode):
        """Constructs a ConnectivityList by wrapping given PropertyField."""
        self._field = field
        self._mode = mode
        self._scoping = scoping

        self.local_scoping = None

    def __getitem__(self, key: int) -> List[int]:
        """Returns a list of indexes or IDs for a given index or ID, see Mode Enum."""
        if self._mode == ReturnMode.IDS:
            return self._get_ids_from_idx(key)
        elif self._mode == ReturnMode.IDX:
            return self._get_idx_from_idx(key)
        raise ValueError(f"ReturnMode has an incorrect value")

    def _get_ids_from_idx(self, idx: int) -> List[int]:
        """Helper method to retrieve list of IDs from a given index."""
        return self._to_ids(self._get_idx_from_idx(idx))

    def _get_idx_from_idx(self, idx: int) -> List[int]:
        """Helper method to retrieve list of indexes from a given index."""
        return self._field.get_entity_data(idx)

    @property
    def by_id(self) -> ConnectivityListById:
        """Returns an equivalent list which accepts an ID instead of an index in __getitem__."""
        return ConnectivityListById(self._field, self._scoping, self._mode)
        
    def _to_ids(self, indices) -> List[int]:
        """Helper method to convert a list of indexes into a list of IDs."""
        if not self.local_scoping:
            self.local_scoping = self._scoping.as_local_scoping()

        to_id = self.local_scoping.id
        return list(map(to_id, indices))
    
    def __contains__(self, l: List[int]) -> bool:
        raise NotImplementedError
    
    def __iter__(self) -> ConnectivityListIterator:
        return ConnectivityListIterator(self)

    def __len__(self) -> int:
        """Returns the number of entities."""
        return len(self._field._get_data_pointer())

class ConnectivityListById(ConnectivityListIdx):
    def __init__(self, field: PropertyField, scoping: Scoping, mode: Mode):
        super().__init__(field, scoping, mode)
    
    def __getitem__(self, key: int) -> List[int]:
        if self._mode == ReturnMode.IDS:
            return self._get_ids_from_id(key)
        elif self._mode == ReturnMode.IDX:
            return self._get_idx_from_id(key)
        raise ValueError(f"ReturnMode has an incorrect value")
    
    def _get_ids_from_id(self, id: int) -> List[int]:
        """Helper method to retrieve list of IDs from a given ID."""
        return self._to_ids(self._get_idx_from_id(id))

    def _get_idx_from_id(self, id: int) -> List[int]:
        """Helper method to retrieve list of indexes from a given ID."""
        return self._field.get_entity_data_by_id(id)
    
    def __contains__(self, l: List[int]) -> bool:
        raise NotImplementedError

    def __iter__(self) -> ConnectivityListIterator:
        return super().__iter__()

    def __len__(self) -> int:
        return len(self._field._get_data_pointer())