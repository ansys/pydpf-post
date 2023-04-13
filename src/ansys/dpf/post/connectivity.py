from __future__ import annotations

from collections.abc import Sequence, Iterator
from enum import Enum
from ansys.dpf.core.property_field import PropertyField
from ansys.dpf.core.scoping import Scoping

from typing import List

class Mode(Enum):
    IDS_FROM_IDX = 1
    IDS_FROM_ID  = 2
    IDX_FROM_IDX = 3
    IDX_FROM_ID  = 4

class ConnectivityList(Sequence):
    def __init__(self, field: PropertyField, scoping: Scoping, mode: Mode):
        self._field = field
        self._mode = mode
        self._scoping = scoping

    def __getitem__(self, key: int) -> List[int]:
        if   self._mode == Mode.IDS_FROM_IDX:
            return self._get_ids_from_idx(key)
        elif self._mode == Mode.IDS_FROM_ID:
            return self._get_ids_from_id(key)
        elif self._mode == Mode.IDX_FROM_IDX:
            return self._get_idx_from_idx(key)
        elif self._mode == Mode.IDX_FROM_ID:
            return self._get_idx_from_id(key)
    
    def _get_ids_from_idx(self, idx: int) -> List[int]:
        return self._to_ids(self._field.get_entity_data(idx))
    
    def _get_ids_from_id(self, id: int) -> List[int]:
        return self._to_ids(self._field.get_entity_data_by_id(id))
    
    def _get_idx_from_idx(self, idx: int) -> List[int]:
        return self._field.get_entity_data(idx)

    def _get_idx_from_id(self, id: int) -> List[int]:
        return self._field.get_entity_data_by_id(id)
    
    @property
    def by_id(self) -> ConnectivityList:
        if self._mode == Mode.IDS_FROM_IDX:
            return ConnectivityList(self._field, self._scoping, Mode.IDS_FROM_ID)
        elif self._mode == Mode.IDX_FROM_IDX:
            return ConnectivityList(self._field, self._scoping, Mode.IDX_FROM_ID)
        return self

    def _to_ids(self, indices) -> List[int]:
        to_id = self._scoping.id
        return list(map(to_id, indices))

    def __len__(self) -> int:
        return len(self._field._get_data_pointer())
        
