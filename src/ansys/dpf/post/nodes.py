"""This module contains Node-related wrapper classes."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from ansys.dpf.core import errors
import ansys.dpf.core.nodes as core_nodes


class Node:
    """Wrapper class around dpf.core.nodes.Node."""

    def __init__(self, node: core_nodes.Node):
        """Constructs a Node from its index and the original list."""
        self._node = node

    @property
    def coordinates(self) -> List[float]:
        """Cartersian coordinates of the node."""
        return self._node.coordinates

    @property
    def id(self) -> int:
        """Returns the ID of the node."""
        return self._node.id

    @property
    def index(self) -> int:
        """Returns the index of the node (zero-based)."""
        return self._node.index

    @property
    def to_element_connectivity(self) -> List[int]:
        """Elements indices connected to the node."""
        return self._node.nodal_connectivity

    def __str__(self) -> str:
        """Returns a string representation of the node."""
        return f"Node(id={self.id}, coordinates={self.coordinates})"

    def __repr__(self) -> str:
        """Returns a string representation of the node."""
        return f"Node(id={self.id})"


class _NodeList(ABC):
    """Iterator class for the NodeList."""

    def __init__(self, nodes: core_nodes.Nodes):
        """Constructs an iterator from a Nodes list."""
        self._nodes = nodes
        self._idx = 0

    def __next__(self) -> Node:
        """Returns the next Node in the list."""
        if self._idx >= len(self._nodes):
            raise StopIteration
        ret = self[self._idx]
        self._idx += 1
        return ret

    def __getitem__(self, index: int) -> Node:
        """Returns a Node at a given index."""
        return Node(self._nodes[index])

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes

    def __repr__(self) -> str:
        """Returns a string representation of a _NodeList object."""
        return f"{self.__class__.__name__}({self}, __len__={len(self)})"

    def _short_list(self) -> str:
        _str = "["
        if self.__len__() > 3:
            _fst = Node(self._nodes[0])
            _lst = Node(self._nodes[len(self) - 1])
            _str += f"{_fst}, ..., {_lst}"
        else:
            node_list = [Node(self._nodes[idx]) for idx in range(len(self))]
            _str += ", ".join(map(repr, node_list))
        _str += "]"
        return _str

    def __str__(self) -> str:
        """Returns a string representation of a _NodeList object."""
        return self._short_list()

    @abstractmethod
    def __iter__(self):  # pragma: no cover
        """Returns an iterator object on the list."""


class NodeListByIndex(_NodeList):
    """Node List object using indexes as input."""

    @property
    def by_id(self) -> NodeListById:
        """Returns an equivalent list which accepts IDs as input."""
        return NodeListById(self._nodes)

    def __iter__(self) -> NodeListByIndex:
        """Returns the object to iterate over."""
        self._idx = 0
        return self

    def __contains__(self, node: Node) -> bool:
        """Checks if the given node is in the list."""
        return len(self) > node.index >= 0


class NodeListById(_NodeList):
    """Node List object using IDs as input."""

    def __getitem__(self, id: int) -> Node:  # pylint: disable=redefined-builtin
        """Access a Node with an ID."""
        idx = self._nodes.scoping.index(id)
        try:
            return super().__getitem__(idx)
        except errors.DPFServerException as e:
            if "node not found" in str(e):
                raise ValueError(f"Node with ID={id} not found in the list.")
            else:
                raise e  # pragma: no cover

    def __contains__(self, node: Node) -> bool:
        """Checks if the given node is in the list."""
        return node.id in self._nodes.scoping.ids

    def __iter__(self) -> NodeListByIndex:
        """Returns the object to iterate over."""
        return NodeListByIndex(self._nodes)
