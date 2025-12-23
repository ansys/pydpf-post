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

import ansys.dpf.core as dpf
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post.common import elemental_properties as elt_prop
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation


@fixture
def meshes(allkindofcomplexity):
    simulation = StaticMechanicalSimulation(allkindofcomplexity)
    return simulation.split_mesh_by_properties(
        properties=[elt_prop.material, elt_prop.element_shape]
    )


def test_meshes_core_object(meshes):
    assert isinstance(meshes._core_object, dpf.MeshesContainer)
    assert (
        elt_prop.element_shape in meshes._core_object.labels
        and elt_prop.material in meshes._core_object.labels
    )


def test_meshes_str(meshes):
    assert str(meshes) == str(meshes._core_object)


def test_meshes_get_item(meshes):
    with pytest.raises(
        ValueError, match="Access to a specific Mesh of a Meshes requires"
    ):
        _ = meshes["test"]
    mesh1 = meshes[1]
    assert isinstance(mesh1, post.Mesh)
    assert len(mesh1.node_ids) == 240
    mesh2 = meshes[{elt_prop.material: 1, elt_prop.element_shape: 0}]
    assert isinstance(mesh2, post.Mesh)
    assert len(mesh2.node_ids) == 240


def test_meshes_plot(meshes):
    _ = meshes.plot()


def test_meshes_select(meshes):
    mesh = meshes.select(mat=1, elshape=0)
    assert isinstance(mesh, post.Mesh)
    assert len(mesh.node_ids) == 240
    mesh_none = meshes.select(mat=22, elshape=42)
    assert mesh_none is None
    meshes_mat = meshes.select(mat=5)
    assert isinstance(meshes_mat, post.Mesh)
    assert len(meshes_mat.node_ids) == 248
    meshes_mat = meshes.select(mat=1)
    # print(meshes_mat._core_object)
    assert isinstance(meshes_mat, post.Meshes)
    assert len(meshes_mat) == 2
