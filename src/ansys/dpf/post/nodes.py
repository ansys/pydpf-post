"""This module contains NodeList class."""

from __future__ import annotations

from typing import List

from collections.abc import Collection, Iterator

import ansys.dpf.core.nodes as nodes

class Node:
    def __init__(self, nodes: nodes.Nodes, index: int):
        self._nodes = nodes
        self._index = index
    
    def _resolve(self) -> nodes.Node:
        return self._nodes[self._index]

    @property
    def index(self) -> int:
        return self._resolve().index
    
    @property
    def id(self) -> int:
        return self._resolve().id
    
    @property
    def coordinates(self) -> List[float]:
        return self._resolve().coordinates
    
    @property
    def nodal_connectivity(self) -> List[int]:
        return self._resolve().nodal_connectivity
    
    def __str__(self) -> str:
        return f"Node(id={self.id}, coordinates={self.coordinates})"
    
    def __repr__(self) -> str:
        return f"Node(id={self.id})"


class NodeListIterator(Iterator):
    """Iterator object for NodeListIdx."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs an iterator from a Node List."""
        self._nodes = nodes
        self._idx = 0

    def __next__(self) -> nodes.Node:
        """Returns the next Node in the list."""
        if self._idx >= self._nodes.__len__():
            raise StopIteration

        ret = Node(self._nodes, self._idx)
        self._idx += 1
        return ret

    def __iter__(self) -> NodeListIterator:
        """Returns a new Iterator object."""
        return NodeListIterator(self._nodes)


class NodeListIdx(Collection):
    """List of Node accessible by index."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs a NodeList from an existing dpf.core.nodes.Nodes object."""
        self._nodes = nodes

    def __getitem__(self, idx: int) -> nodes.Node:
        """Returns a Node at a given index."""
        return Node(self._nodes, idx)

    def __contains__(self, node: nodes.Node) -> bool:
        """Checks if given node is in the list."""
        return node.index >= 0 and node.index < self.__len__()

    def __iter__(self) -> NodeListIterator:
        """Returns an iterator object on the list."""
        return NodeListIterator(self._nodes)

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes

    @property
    def by_id(self) -> NodeListById:
        """Returns an equivalent list Accessible by ID."""
        return NodeListById(self._nodes)
    
    def _short_list(self) -> str:
        _str = "["
        if self.__len__() > 3:
            _fst = Node(self._nodes, 0)
            _lst = Node(self._nodes, self.__len__()-1)
            _str += f"{_fst}, ..., {_lst}"
        else:
            el_list = [Node(self._nodes, idx) for idx in range(self.__len__())]
            _str += ", ".join(map(lambda el: repr(el), el_list))
        _str += "]"
        return _str
    
    def __str__(self) -> str:
        return self._short_list()
    
    def __repr__(self) -> str:
        return f"NodeListIdx({self.__str__()}, __len__={self.__len__()})"


class NodeListById(NodeListIdx):
    """List of node accessible by ID."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs a list from an existing core.nodes.Nodes object."""
        super().__init__(nodes)

    def __getitem__(self, id: int) -> nodes.Node:
        """Returns a Node for a given ID."""
        idx = self._nodes.scoping.index(id)
        return super().__getitem__(idx)

    def __contains__(self, node: nodes.Node) -> bool:
        """Checks if the given node is in the list."""
        return node.id in self._nodes.scoping.ids

    def __iter__(self) -> NodeListIterator:
        """Returns an iterator object on the list."""
        return super().__iter__()

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes

    def __repr__(self) -> str:
        return f"NodeListById({super().__str__()}, __len__={self.__len__()})"