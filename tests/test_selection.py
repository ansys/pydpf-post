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

import numpy as np
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.selection import SpatialSelection
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1,
)


def test_spatial_selection_select_nodes(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection._selection.progress_bar = False
    selection.select_nodes([1, 2, 3])
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.locations.nodal
    assert np.allclose(scoping.ids, [1, 2, 3])


def test_spatial_selection_select_elements(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection._selection.progress_bar = False
    selection.select_elements([1, 2, 3, 4])
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.locations.elemental
    assert np.allclose(scoping.ids, [1, 2, 3, 4])
    ids = selection.apply_to(simulation)
    assert np.allclose(ids, [1, 2, 3, 4])


def test_spatial_selection_select_named_selection(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection._selection.progress_bar = False
    selection.select_named_selection(
        simulation.mesh.named_selections.keys()[0],
        location=post.selection.locations.nodal,
    )
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.locations.nodal
    assert scoping.ids.size == 12970
    assert 1857 in scoping.ids
    assert 14826 in scoping.ids
    ids = selection.apply_to(simulation)
    assert len(ids) == 12970
    assert 1857 in ids
    assert 14826 in ids


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Faces added with ansys-dpf-server 2024.1.pre0.",
)
class TestSpatialSelectionFaces:
    @fixture
    def fluent_simulation(self):
        fluid_example_files = examples.download_fluent_axial_comp()
        ds = dpf.DataSources()
        ds.set_result_file_path(
            fluid_example_files["cas"][0],
            key="cas",
        )
        ds.add_file_path(
            fluid_example_files["dat"][0],
            key="dat",
        )
        return post.FluidSimulation(ds)  # noqa

    def test_spatial_selection_select_faces(self, fluent_simulation):
        selection = SpatialSelection()
        selection._selection.progress_bar = False
        selection.select_faces(fluent_simulation.mesh.face_ids)
        scoping = selection._evaluate_on(fluent_simulation)
        assert scoping.location == post.locations.faces
        assert np.allclose(scoping.ids, fluent_simulation.mesh.face_ids)

    def test_spatial_selection_select_nodes_of_faces(self, fluent_simulation):
        selection = SpatialSelection()
        selection._selection.progress_bar = False
        face_0 = fluent_simulation.mesh.faces[0]
        selection.select_nodes_of_faces(
            faces=[face_0.id],
            mesh=fluent_simulation.mesh,
        )
        scoping = selection._evaluate_on(fluent_simulation)
        assert scoping.location == post.locations.nodal
        assert np.allclose(scoping.ids, face_0.node_ids)

    def test_spatial_selection_select_faces_of_elements(self, fluent_simulation):
        selection = SpatialSelection()
        selection._selection.progress_bar = False
        elem_0 = fluent_simulation.mesh.elements[0]
        selection.select_faces_of_elements(
            elements=[elem_0.id],
            mesh=fluent_simulation.mesh,
        )
        scoping = selection._evaluate_on(fluent_simulation)
        assert scoping.location == post.locations.faces
        if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1:
            list_ref = [11479, 11500, -1, 11502, 11503]
        else:
            list_ref = [12481, 12502, 39941, 43681, 12504, 12505]
        assert np.allclose(scoping.ids, list_ref)


#
#
# def test_spatial_selection_intersect(allkindofcomplexity):
#     solution = post.load_solution(allkindofcomplexity, legacy=False)
#     selection1 = SpatialSelection()
#     selection1.select_nodes([1, 2, 3])
#     _ = selection1._evaluate_on(solution)
#
#     selection = SpatialSelection()
#     selection.select_nodes([1, 2, 3, 4])
#     selection.intersect(selection1)
#     scoping = selection._evaluate_on(solution)
#     assert scoping.location == post.selection.locations.nodal
#     assert np.allclose(scoping.ids, [1, 2, 3])
#     ids = selection.apply_to(solution)
#     assert np.allclose(ids, [1, 2, 3])
