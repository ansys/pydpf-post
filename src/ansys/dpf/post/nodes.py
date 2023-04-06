"""This module contains NodeList class."""

from __future__ import annotations

from collections.abc import Sequence

import ansys.dpf.core.nodes as nodes


class NodeListIdx(Sequence):
    """List of Node."""

    def __init__(self, nodes: nodes.Nodes):
        """Constructs a NodeList from an existing dpf.core.nodes.Nodes object."""
        self._nodes = nodes

    def __getitem__(self, idx: int) -> nodes.Node:
        """Delegates to node_by_id() if by_id, otherwise to node_by_index()."""
        return self._nodes.node_by_index(idx)

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes

    @property
    def by_id(self) -> NodeListById:
        return NodeListById(self._nodes)
    
class NodeListById(NodeListIdx):
    def __init__(self, nodes: nodes.Nodes):
        super().__init__(nodes)

    def __getitem__(self, id: int) -> nodes.Node:
        return self._nodes.node_by_id(id)