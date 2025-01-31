# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
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
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf.post import StaticMechanicalSimulation
from ansys.dpf.post.named_selection import NamedSelection
import conftest


@fixture
def mesh(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.mesh


def test_named_selections(mesh):
    ns = mesh.named_selections
    assert len(ns) == 1
    assert len(ns.keys()) == 1
    assert all([isinstance(n, str) for n in ns])
    assert ns.keys()[0] == "_FIXEDSU"
    n = ns[ns.keys()[0]]
    assert isinstance(n, NamedSelection)
    assert len(n.ids) == 21
    with pytest.raises(KeyError, match="could not be found"):
        _ = ns["test"]


def test_named_selection():
    scoping = dpf.Scoping(ids=[1, 2, 3], location=dpf.locations.nodal)
    ns = NamedSelection(name="test", scoping=scoping)
    assert ns._scoping == scoping
    assert ns.name == "test"
    assert ns.id(0) == 1
    assert ns.index(1) == 0
    assert ns.location == dpf.locations.nodal
    assert ns.size == 3
    ref = """NamedSelection 'test'
 with DPF  Scoping: 
  with Nodal location and 3 entities
"""  # noqa
    assert repr(ns) == ref

    with pytest.raises(ValueError, match="No list of entity IDs given."):
        _ = NamedSelection(name="test")
    with pytest.raises(ValueError, match="NamedSelection accepts only"):
        _ = NamedSelection(name="test", node_ids=[0, 1], element_ids=[0, 1])

    node_ns = NamedSelection(name="test", node_ids=[0, 1])
    assert node_ns._scoping.location == dpf.locations.nodal
    element_ns = NamedSelection(name="test", element_ids=[0, 1])
    assert element_ns._scoping.location == dpf.locations.elemental
    if conftest.SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0:
        face_ns = NamedSelection(name="test", face_ids=[0, 1])
        assert face_ns._scoping.location == dpf.locations.faces
    cell_ns = NamedSelection(name="test", cell_ids=[0, 1])
    assert cell_ns._scoping.location == dpf.locations.elemental
