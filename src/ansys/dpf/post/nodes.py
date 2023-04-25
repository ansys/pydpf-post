"""This module contains NodeList class."""

from __future__ import annotations

from collections.abc import Collection, Iterator

import ansys.dpf.core.nodes as nodes

class NodeListIterator(Iterator):
    def __init__(self, nodes: NodeListIdx):
        self._nodes = nodes
        self._idx = 0

    def __next__(self) -> nodes.Node:
        if self._idx >= self._nodes.__len__():
            raise StopIteration	
        
        ret = self._nodes[self._idx]
        self._idx += 1
        return ret

    def __iter__(self) -> NodeListIterator:
        return NodeListIterator(self._nodes)

class NodeListIdx(Collection):
    """List of Node accessible by index."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs a NodeList from an existing dpf.core.nodes.Nodes object."""
        self._nodes = nodes

    def __getitem__(self, idx: int) -> nodes.Node:
        """Returns a Node at a given index."""
        return self._nodes.node_by_index(idx)

    def __contains__(self, node: nodes.Node) -> bool:
        return node.index >= 0 and node.index < self.__len__()

    def __iter__(self) -> NodeListIterator:
        return NodeListIterator(self)

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes

    @property
    def by_id(self) -> NodeListById:
        """Returns an equivalent list Accessible by ID."""
        return NodeListById(self._nodes)


class NodeListById(NodeListIdx):
    """List of node accessible by ID."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs a list from an existing core.nodes.Nodes object."""
        super().__init__(nodes)

    def __getitem__(self, id: int) -> nodes.Node:
        """Returns a Node for a given ID."""
        return self._nodes.node_by_id(id)

    def __contains__(self, node: nodes.Node) -> bool:
        return node.id in self._nodes.scoping.ids
    
    def __iter__(self) -> NodeListIterator:
        return super().__iter__()

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes
