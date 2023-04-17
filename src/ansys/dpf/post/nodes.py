"""This module contains NodeList class."""

from __future__ import annotations

from collections.abc import Sequence

import ansys.dpf.core.nodes as nodes


class NodeListIdx(Sequence):
    """List of Node accessible by index."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs a NodeList from an existing dpf.core.nodes.Nodes object."""
        self._nodes = nodes

    def __getitem__(self, idx: int) -> nodes.Node:
        """Returns a Node at a given index."""
        return self._nodes.node_by_index(idx)

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

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes
