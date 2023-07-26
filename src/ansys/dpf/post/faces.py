"""This module contains Face and Faces classes."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from ansys.dpf.core import errors
from ansys.dpf.core import faces as core_faces

from ansys.dpf.post.elements import ElementType
from ansys.dpf.post.nodes import NodeListByIndex


class Face:
    """Proxy class wrapping dpf.core.faces.Face."""

    def __init__(self, face: core_faces.Face):
        """Constructs a Proxy Face object."""
        self._face = face

    @property
    def node_ids(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.faces.Face.node_ids`."""
        return self._face.node_ids

    @property
    def id(self) -> int:
        """See :py:meth:`ansys.dpf.core.faces.Face.id`."""
        return self._face.id

    @property
    def index(self) -> int:
        """See :py:meth:`ansys.dpf.core.faces.Face.index`."""
        return self._face.index

    @property
    def nodes(self) -> NodeListByIndex:
        """See :py:meth:`ansys.dpf.post.nodes.NodeListByIndex`."""
        return NodeListByIndex(self._face.nodes)

    @property
    def num_nodes(self) -> int:
        """See :py:meth:`ansys.dpf.core.faces.Face.n_nodes`."""
        return self._face.n_nodes

    @property
    def type_info(self) -> ElementType:
        """Gets an element descriptor, See :py:meth:`ansys.dpf.core.faces.Face.id`."""
        return ElementType(self._face.type.value)

    @property
    def type(self) -> core_elements.element_types:
        """Returns the Element Type."""
        return self._face.type

    @property
    def to_node_connectivity(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.faces.Face.connectivity`."""
        return self._face.connectivity

    def __repr__(self) -> str:
        """Returns string representation of a Face."""
        return f"Face(type={self.type},index={self.index},id={self.id})"

    def __str__(self) -> str:
        """Returns string representation of a Face."""
        return str(self._face)


class _FaceList(ABC):
    """Iterator class for the FaceList."""

    def __init__(self, face_list: core_faces.Faces):
        """Constructs an Iterator from a face list."""
        self._face_list = face_list
        self._idx = 0

    def __next__(self) -> Face:
        """Returns the next Face in the list."""
        if self._idx >= len(self._face_list):
            raise StopIteration
        ret = self[self._idx]
        self._idx += 1
        return ret

    def __getitem__(self, index: int) -> Face:
        """Returns a post.Face based on an index in the current list."""
        return Face(self._face_list[index])

    def __len__(self) -> int:
        """Returns the number of faces in the list."""
        return self._face_list.n_faces

    def __repr__(self) -> str:
        """Returns a string representation of _FaceList object."""
        return f"{self.__class__.__name__}({self}, __len__={len(self)})"

    def _short_list(self) -> str:
        _str = "["
        if self.__len__() > 3:
            _fst = Face(self._face_list[0]).type_info.name
            _lst = Face(self._face_list[len(self) - 1]).type_info.name
            _str += f"{_fst}, ..., {_lst}"
        else:
            face_list = [Face(self._face_list[idx]) for idx in range(len(self))]
            _str += ", ".join(map(lambda el: el.type_info.name, face_list))
        _str += "]"
        return _str

    def __str__(self) -> str:
        """Returns a string representation of a _FaceList object."""
        return self._short_list()

    @abstractmethod
    def __iter__(self):  # pragma: no cover
        """Returns the object to iterate on."""
        raise NotImplementedError


class FaceListByIndex(_FaceList):
    """Face list object using indexes as input."""

    @property
    def by_id(self) -> FaceListById:
        """Returns an equivalent list which accepts IDs as input."""
        return FaceListById(self._face_list)

    def __iter__(self) -> FaceListByIndex:
        """Returns the object to iterate over."""
        self._idx = 0
        return self

    def __contains__(self, face: Face) -> bool:
        """Checks if the given element in the list."""
        return len(self) > face.index >= 0


class FaceListById(_FaceList):
    """Face list object using IDs as input."""

    def __getitem__(self, id: int) -> Face:  # pylint: disable=redefined-builtin
        """Access a Face with an ID."""
        idx = self._face_list.scoping.index(id)
        try:
            return super().__getitem__(idx)
        except errors.DPFServerException as e:
            if "face not found" in str(e):
                raise ValueError(f"Face with ID={id} not found in the list.")
            else:
                raise e  # pragma: no cover

    def __contains__(self, face: Face) -> bool:
        """Checks if the given face is in the list."""
        return face.id in self._face_list.scoping.ids

    def __iter__(self) -> FaceListByIndex:
        """Returns the object to iterate over."""
        return FaceListByIndex(self._face_list)
