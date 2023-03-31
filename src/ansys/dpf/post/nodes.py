"""This module contains NodeList class."""

from __future__ import annotations

from collections.abc import Sequence

import ansys.dpf.core.nodes as nodes


class NodeList(Sequence):
    """List of Node."""

    def __init__(self, nodes: nodes.Nodes, by_id=True):
        """Constructs a NodeList from an existing dpf.core.nodes.Nodes object."""
        self._nodes = nodes
        self.by_id = by_id

    def __getitem__(self, key: int) -> nodes.Node:
        """Delegates to node_by_id() if by_id, otherwise to node_by_index()."""
        if self.by_id:
            return self._nodes.node_by_id(key)
        else:
            return self._nodes.node_by_index(key)

    def __len__(self) -> int:
        """Returns the number of nodes in the list."""
        return self._nodes.n_nodes
