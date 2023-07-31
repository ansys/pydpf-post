"""Module containing wrapper class for the connectivities property fields."""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import List

from ansys.dpf.core.property_field import PropertyField
from ansys.dpf.core.scoping import Scoping


class ReturnMode(Enum):
    """Enum made for internal use, to dictate the behavior of _ConnectivityList."""

    IDS = 1
    IDX = 2


class _ConnectivityList(ABC):
    """Abstract class for ConnectivityList objects."""

    def __init__(
        self,
        field: PropertyField,
        mode: ReturnMode,
        scoping: Scoping,
    ):
        """Constructs a _ConnectivityList by wrapping the given PropertyField.

        Parameters
        ----------
        field:
            Field of connectivity.
        mode:
            Whether to return indexes or IDs.
        scoping:
            Element or node scoping to map returned indexes to IDs.
        """
        self._field = field
        if mode not in [ReturnMode.IDS, ReturnMode.IDX]:
            raise ValueError("'mode' argument must be a valid ReturnMode value")
        self._mode = mode
        self._scoping = scoping
        self._idx = 0
        self.local_scoping = None

    def __next__(self) -> List[int]:
        """Returns the next element in the list."""
        if self._idx >= len(self):
            raise StopIteration
        out = self.__getitem__(self._idx)
        self._idx += 1
        return out

    def __getitem__(self, idx: int) -> List[int]:
        """Returns, for a given entity index, the connected indexes or IDs (see ReturnMode)."""
        if self._mode == ReturnMode.IDS:
            return self._get_ids_from_idx(idx)
        elif self._mode == ReturnMode.IDX:
            return self._get_idx_from_idx(idx)

    def __len__(self) -> int:
        """Returns the number of entities."""
        return self._field.scoping.size

    def __repr__(self) -> str:
        """Returns string representation of a _ConnectivityList object."""
        return f"{self.__class__.__name__}({self.__str__()}, __len__={self.__len__()})"

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
        """Returns string representation of a _ConnectivityList object."""
        return self._short_list()

    def _get_ids_from_idx(self, idx: int) -> List[int]:
        """Helper method to retrieve list of IDs from a given index."""
        return self._to_ids(self._get_idx_from_idx(idx))

    def _get_idx_from_idx(self, idx: int) -> List[int]:
        """Helper method to retrieve list of indexes from a given index."""
        return self._field.get_entity_data(idx).tolist()

    def _to_ids(self, indices) -> List[int]:
        """Helper method to convert a list of indexes into a list of IDs."""
        if not self.local_scoping:
            self.local_scoping = self._scoping.as_local_scoping()

        to_id = self.local_scoping.id
        return list(map(to_id, indices))

    @abstractmethod
    def __iter__(self):  # pragma: no cover
        """Returns the object to iterate on."""


class ConnectivityListByIndex(_ConnectivityList):
    """Connectivity list object using indexes as input."""

    @property
    def by_id(self) -> ConnectivityListById:
        """Returns an equivalent list which accepts IDs as input."""
        return ConnectivityListById(
            field=self._field, mode=self._mode, scoping=self._scoping
        )

    def __iter__(self) -> ConnectivityListByIndex:
        """Returns the object to iterate on."""
        self._idx = 0
        return self


class ConnectivityListById(_ConnectivityList):
    """Connectivity list object using IDs as input."""

    def __getitem__(self, id: int) -> List[int]:  # pylint: disable=redefined-builtin
        """Returns, for a given entity ID, the connected indexes or IDs (see ReturnMode)."""
        idx = self._field.scoping.index(id)
        return super().__getitem__(idx)

    def __iter__(self) -> ConnectivityListByIndex:
        """Returns the object to iterate on."""
        return ConnectivityListByIndex(
            field=self._field,
            mode=self._mode,
            scoping=self._scoping,
        )
