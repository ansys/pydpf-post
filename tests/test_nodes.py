# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import pytest

from ansys.dpf import core as dpf
from ansys.dpf.post import examples
from ansys.dpf.post.nodes import Node, NodeListById, NodeListByIndex


def test_node():
    model = dpf.Model(examples.find_static_rst())
    core_nodes = model.metadata.meshed_region.nodes
    node = Node(node=core_nodes[0])
    ref = [0.015, 0.045, 0.015]
    assert node.coordinates == ref
    assert node.id == 1
    assert node.index == 0
    ref = [0, 1, 2, 3, 4, 5, 6, 7]
    assert list(node.to_element_connectivity) == ref
    ref = "Node(id=1, coordinates=[0.015, 0.045, 0.015])"
    assert str(node) == ref
    ref = "Node(id=1)"
    assert repr(node) == ref


def test_nodes_nodes_list_by_idx():
    model = dpf.Model(examples.find_static_rst())
    core_nodes = model.metadata.meshed_region.nodes
    nodes_list_by_index = NodeListByIndex(nodes=core_nodes)
    for i in nodes_list_by_index:
        assert isinstance(i, Node)
    for i in nodes_list_by_index:
        assert isinstance(i, Node)
    assert nodes_list_by_index[1].id == 2
    assert len(nodes_list_by_index) == 81
    ref = (
        "NodeListByIndex([Node(id=1, coordinates=[0.015, 0.045, 0.015]), ..., "
        "Node(id=81, coordinates=[0.03, 0.045, 0.0075])], __len__=81)"
    )
    assert repr(nodes_list_by_index) == ref
    ref = (
        "[Node(id=1, coordinates=[0.015, 0.045, 0.015]), ..., "
        "Node(id=81, coordinates=[0.03, 0.045, 0.0075])]"
    )
    assert str(nodes_list_by_index) == ref
    assert nodes_list_by_index[0] in nodes_list_by_index
    nodes_list_by_id = nodes_list_by_index.by_id
    assert isinstance(nodes_list_by_id, NodeListById)


def test_nodes_nodes_list_by_id():
    model = dpf.Model(examples.find_static_rst())
    core_nodes = model.metadata.meshed_region.nodes
    nodes_list_by_id = NodeListById(nodes=core_nodes)
    for i in nodes_list_by_id:
        assert isinstance(i, Node)
    assert isinstance(nodes_list_by_id[5], Node)
    assert nodes_list_by_id[5] in nodes_list_by_id
    with pytest.raises(ValueError, match="not found"):
        _ = nodes_list_by_id[0]
