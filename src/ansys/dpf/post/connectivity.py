"""Module containing wrapper class for the connectivities property fields."""

from __future__ import annotations

from collections.abc import Collection, Iterator
from enum import Enum
from typing import List

from ansys.dpf.core.any import Any
from ansys.dpf.core.property_field import PropertyField
from ansys.dpf.core.scoping import Scoping


class ReturnMode(Enum):
    """Enum made for internal use, to dictate the behavior of ConnectivityList."""

    IDS = 1
    IDX = 2


class ConnectivityListIterator(Iterator):
    """Iterator class for Connectivity Lists."""

    def __init__(self, conn_list: PropertyField):
        """Constructs an iterator from an existing list."""
        self._conn_list = conn_list
        self._idx = 0

    def __next__(self) -> Any:
        """Returns the next element in the list."""
        if self._idx >= self._conn_list.__len__():
            raise StopIteration

        ret = self._conn_list.get_entity_data(self._idx)
        self._idx += 1
        return ret

    def __iter__(self) -> ConnectivityListIterator:
        """Returns a new ConnectivityListIterator."""
        return ConnectivityListIdx(field=self._conn_list)


class ConnectivityListIdx(Collection):
    """Very basic wrapper around elemental and nodal connectivities fields."""

    def __init__(
        self,
        field: PropertyField,
        scoping: Union[Scoping, None] = None,
        mode: Union[ReturnMode, None] = None,
    ):
        """Constructs a ConnectivityList by wrapping given PropertyField."""
        self._field = field
        self._mode = mode
        self._scoping = scoping
        self._idx = 0
        self.local_scoping = None

    def __next__(self) -> Any:
        """Returns the next element in the list."""
        if self._idx >= len(self):
            raise StopIteration
        if self._mode == ReturnMode.IDS:
            out = self._get_ids_from_idx(self._idx)
        elif self._mode == ReturnMode.IDX:
            out = self._get_idx_from_idx(self._idx)
        else:
            raise ValueError(f"ReturnMode has an incorrect value")
        self._idx += 1
        return out

    def __getitem__(self, idx: int) -> List[int]:
        """Returns a list of indexes or IDs for a given index, see ReturnMode Enum."""
        if self._mode == ReturnMode.IDS:
            return self._get_ids_from_idx(idx)
        elif self._mode == ReturnMode.IDX:
            return self._get_idx_from_idx(idx)
        raise ValueError("ReturnMode has an incorrect value")

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
        """Not implemented."""
        raise NotImplementedError

    def __iter__(self) -> ConnectivityListIterator:
        """Returns an iterator object on the list."""
        return ConnectivityListIterator(self._field)

    def __len__(self) -> int:
        """Returns the number of entities."""
        return len(self._field._get_data_pointer())

    def _short_list(self) -> str:
        _str = "["
        if self.__len__() > 3:
            _fst = self._field.get_entity_data(0)
            _lst = self._field.get_entity_data(self.__len__() - 1)
            if self._mode == ReturnMode.IDS:
                _fst = self._to_ids(_fst)
                _lst = self._to_ids(_lst)
            _str += f"{_fst}, ..., {_lst}"
        else:
            conn_list = [
                self._field.get_entity_data(idx) for idx in range(self.__len__())
            ]
            if self._mode == ReturnMode.IDS:
                conn_list = list(map(self._to_ids, conn_list))
            _str += ", ".join(map(str, conn_list))
        _str += "]"
        return _str

    def __str__(self) -> str:
        """Returns string representation of ConnectivityListIdx."""
        return self._short_list()

    def __repr__(self) -> str:
        """Returns string representation of ConnectivityListIdx."""
        return f"ConnectivityListIdx({self.__str__()},__len__={self.__len__()})"


class ConnectivityListById(ConnectivityListIdx):
    """Connectivity List indexed by ID."""

    def __init__(self, field: PropertyField, scoping: Scoping, mode: ReturnMode):
        """Constructs a Connectivity list from a given PropertyField."""
        super().__init__(field, scoping, mode)

    def __getitem__(self, id: int) -> List[int]:  # noqa: W0622
        """Returns a list of indexes or IDs for a given ID, see ReturnMode Enum."""
        idx = self._field.scoping.index(id)
        return super().__getitem__(idx)

    def __contains__(self, l: List[int]) -> bool:
        """Not Implemented."""
        raise NotImplementedError

    def __iter__(self) -> ConnectivityListIterator:
        """Returns an iterator object on the list."""
        return super().__iter__()

    def __len__(self) -> int:
        """Returns the number of entities."""
        return len(self._field._get_data_pointer())

    def __repr__(self) -> str:
        """String representation of ConnectivityListById."""
        return f"ConnectivityListById({super().__str__()}, __len__={self.__len__()})"
