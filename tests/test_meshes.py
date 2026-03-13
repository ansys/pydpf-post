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

import ansys.dpf.core as dpf
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post.common import elemental_properties as elt_prop
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0


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
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        mesh2 = meshes[{elt_prop.material: 1, elt_prop.element_shape: 1}]
    else:
        mesh2 = meshes[{elt_prop.material: 1, elt_prop.element_shape: 0}]
    assert isinstance(mesh2, post.Mesh)
    assert len(mesh2.node_ids) == 240


def test_meshes_plot(meshes):
    _ = meshes.plot()


def test_meshes_select(meshes):
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        mesh = meshes.select(mat=1, elshape=1)
    else:
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


def test_meshes_select_solids_all(meshes):
    """Test select_solids with no additional filters returns all solid meshes."""
    result = meshes.select_solids()
    assert result is not None
    assert isinstance(result, post.Meshes)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=2 (SOLID): meshes 2 (mat=1) and 6 (mat=2) → 2 meshes
        assert len(result) == 2
    else:
        # elshape=1 (SOLID): meshes 2 (mat=1) and 8 (mat=2) → 2 meshes
        assert len(result) == 2


def test_meshes_select_solids_with_material_filter_single(meshes):
    """Test select_solids filtered by mat=1 returns a single solid mesh."""
    result = meshes.select_solids(mat=1)
    assert result is not None
    assert isinstance(result, post.Mesh)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=2, mat=1 → mesh 2: 1856 nodes, 343 elements
        assert len(result.node_ids) == 1856
        assert len(result.element_ids) == 343
    else:
        # elshape=1, mat=1 → mesh 2: 1856 nodes, 343 elements
        assert len(result.node_ids) == 1856
        assert len(result.element_ids) == 343


def test_meshes_select_solids_with_material_filter_mat2(meshes):
    """Test select_solids filtered by mat=2 returns a single solid mesh."""
    result = meshes.select_solids(mat=2)
    assert result is not None
    assert isinstance(result, post.Mesh)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=2, mat=2 → mesh 6: 12970 nodes, 8709 elements
        assert len(result.node_ids) == 12970
        assert len(result.element_ids) == 8709
    else:
        # elshape=1, mat=2 → mesh 8: 12970 nodes, 8709 elements
        assert len(result.node_ids) == 12970
        assert len(result.element_ids) == 8709


def test_meshes_select_solids_with_material_filter_multiple(meshes):
    """Test select_solids filtered by a list of materials."""
    # Both mat=1 and mat=2 have solids in both server versions
    result = meshes.select_solids(mat=[1, 2])
    assert result is not None
    assert isinstance(result, post.Meshes)
    assert len(result) == 2


def test_meshes_select_solids_nonexistent_material(meshes):
    """Test select_solids with a material that does not exist returns None."""
    result = meshes.select_solids(mat=9999)
    assert result is None


def test_meshes_select_solids_ignores_unknown_kwargs(meshes):
    """Test that unknown kwargs are ignored and do not cause errors."""
    result_all = meshes.select_solids()
    result_with_unknown = meshes.select_solids(nonexistent_label=42)
    assert isinstance(result_all, post.Meshes)
    assert isinstance(result_with_unknown, post.Meshes)
    assert len(result_with_unknown) == len(result_all)


def test_meshes_select_solids_with_list_material_partial_match(meshes):
    """Test select_solids with a list of materials where only some match."""
    # mat=1 has solids, mat=9999 does not → 1 mesh
    result = meshes.select_solids(mat=[1, 9999])
    assert result is not None
    assert isinstance(result, post.Mesh)
    assert len(result.node_ids) == 1856


def test_meshes_select_solids_no_match_for_material_without_solids(meshes):
    """Test select_solids with a material that exists but has no solid elements."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # mat=5 has elshape=1 (SHELL), no elshape=2 (SOLID)
        result = meshes.select_solids(mat=5)
        assert result is None
    else:
        # mat=5 has elshape=0 (SHELL), no elshape=1 (SOLID)
        result = meshes.select_solids(mat=5)
        assert result is None


def test_meshes_select_shells_all(meshes):
    """Test select_shells with no additional filters returns all shell meshes."""
    result = meshes.select_shells()
    assert result is not None
    assert isinstance(result, post.Meshes)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=1 (SHELL): meshes 0,1,3,5,13 → 5 meshes
        assert len(result) == 5
    else:
        # elshape=0 (SHELL): meshes 0,1,3,6,9,15,17 → 7 meshes
        assert len(result) == 7


def test_meshes_select_shells_with_material_filter_single(meshes):
    """Test select_shells filtered by a specific material returns a single mesh."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=1, mat=5 → mesh 0: 248 nodes, 105 elements
        result = meshes.select_shells(mat=5)
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 248
        assert len(result.element_ids) == 105
    else:
        # elshape=0, mat=5 → mesh 0: 248 nodes, 105 elements
        result = meshes.select_shells(mat=5)
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 248
        assert len(result.element_ids) == 105


def test_meshes_select_shells_with_material_filter_mat1(meshes):
    """Test select_shells filtered by mat=1."""
    result = meshes.select_shells(mat=1)
    assert result is not None
    assert isinstance(result, post.Mesh)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=1, mat=1 → mesh 1: 240 nodes, 203 elements
        assert len(result.node_ids) == 240
        assert len(result.element_ids) == 203
    else:
        # elshape=0, mat=1 → mesh 1: 240 nodes, 203 elements
        assert len(result.node_ids) == 240
        assert len(result.element_ids) == 203


def test_meshes_select_shells_with_material_filter_multiple(meshes):
    """Test select_shells filtered by a list of materials."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=1, mat=[1,5] → mesh 1 (mat=1) and mesh 0 (mat=5) → 2 meshes
        result = meshes.select_shells(mat=[1, 5])
        assert result is not None
        assert isinstance(result, post.Meshes)
        assert len(result) == 2
    else:
        # elshape=0, mat=[6,8] → mesh 9 (mat=6) and mesh 17 (mat=8) → 2 meshes
        result = meshes.select_shells(mat=[6, 8])
        assert result is not None
        assert isinstance(result, post.Meshes)
        assert len(result) == 2


def test_meshes_select_shells_nonexistent_material(meshes):
    """Test select_shells with a material that does not exist returns None."""
    result = meshes.select_shells(mat=9999)
    assert result is None


def test_meshes_select_shells_ignores_unknown_kwargs(meshes):
    """Test that unknown kwargs are ignored and do not cause errors."""
    result_all = meshes.select_shells()
    result_with_unknown = meshes.select_shells(nonexistent_label=42)
    assert isinstance(result_all, post.Meshes)
    assert isinstance(result_with_unknown, post.Meshes)
    assert len(result_with_unknown) == len(result_all)


def test_meshes_select_shells_no_match_for_material_without_shells(meshes):
    """Test select_shells with a material that exists but has no shell elements."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # mat=2 has elshape=2 (SOLID), no elshape=1 (SHELL)
        result = meshes.select_shells(mat=2)
        assert result is None
    else:
        # mat=2 has elshape=1 (SOLID), no elshape=0 (SHELL)
        result = meshes.select_shells(mat=2)
        assert result is None


def test_meshes_select_shells_specific_node_and_element_counts(meshes):
    """Test that select_shells returns meshes with expected node/element counts."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # mat=7, elshape=1 → mesh 13: 553 nodes, 276 elements
        result = meshes.select_shells(mat=7)
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 553
        assert len(result.element_ids) == 276
    else:
        # mat=7, elshape=0 → mesh 15: 553 nodes, 276 elements
        result = meshes.select_shells(mat=7)
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 553
        assert len(result.element_ids) == 276


def test_meshes_select_beams_all(meshes):
    """Test select_beams with no additional filters returns all beam meshes."""
    result = meshes.select_beams()
    assert result is not None
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=3 (BEAM): mesh 4 (mat=4) → 1 mesh (single Mesh)
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 17
        assert len(result.element_ids) == 16
    else:
        # elshape=2 (BEAM): meshes 5,10,12,14,16 → 5 meshes
        assert isinstance(result, post.Meshes)
        assert len(result) == 5


def test_meshes_select_beams_with_material_filter_single(meshes):
    """Test select_beams filtered by a specific material returns a single mesh."""
    result = meshes.select_beams(mat=4)
    assert result is not None
    assert isinstance(result, post.Mesh)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=3, mat=4 → mesh 4: 17 nodes, 16 elements
        assert len(result.node_ids) == 17
        assert len(result.element_ids) == 16
    else:
        # elshape=2, mat=4 → mesh 5: 17 nodes, 16 elements
        assert len(result.node_ids) == 17
        assert len(result.element_ids) == 16


def test_meshes_select_beams_nonexistent_material(meshes):
    """Test select_beams with a material that does not exist returns None."""
    result = meshes.select_beams(mat=9999)
    assert result is None


def test_meshes_select_beams_ignores_unknown_kwargs(meshes):
    """Test that unknown kwargs are ignored and do not cause errors."""
    result_all = meshes.select_beams()
    result_with_unknown = meshes.select_beams(nonexistent_label=42)
    if isinstance(result_all, post.Mesh):
        assert isinstance(result_with_unknown, post.Mesh)
        assert len(result_all.node_ids) == len(result_with_unknown.node_ids)
    else:
        assert isinstance(result_with_unknown, post.Meshes)
        assert len(result_with_unknown) == len(result_all)


def test_meshes_select_beams_no_match_for_material_without_beams(meshes):
    """Test select_beams with a material that exists but has no beam elements."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # mat=1 has elshape=1 (SHELL) and elshape=2 (SOLID), no elshape=3 (BEAM)
        result = meshes.select_beams(mat=1)
        assert result is None
    else:
        # mat=1 has elshape=0 (SHELL) and elshape=1 (SOLID), no elshape=2 (BEAM)
        result = meshes.select_beams(mat=1)
        assert result is None


def test_meshes_select_beams_with_list_material(meshes):
    """Test select_beams with a list of material values."""
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=3: only mat=4 exists → 1 mesh
        result = meshes.select_beams(mat=[4, 9999])
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 17
    else:
        # elshape=2: mat=4 and mat=12 → 2 meshes
        result = meshes.select_beams(mat=[4, 12])
        assert result is not None
        assert isinstance(result, post.Meshes)
        assert len(result) == 2


def test_meshes_select_beams_legacy_specific_materials(meshes):
    """Test select_beams with specific legacy materials that have beam elements."""
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # elshape=2, mat=12 → mesh 10: 2 nodes, 1 element
        result = meshes.select_beams(mat=12)
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 2
        assert len(result.element_ids) == 1

        # elshape=2, mat=15 → mesh 12: 2 nodes, 1 element
        result = meshes.select_beams(mat=15)
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 2
        assert len(result.element_ids) == 1

        # elshape=2, mat=16 → mesh 14: 2 nodes, 1 element
        result = meshes.select_beams(mat=16)
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 2
        assert len(result.element_ids) == 1

        # elshape=2, mat=17 → mesh 16: 2 nodes, 1 element
        result = meshes.select_beams(mat=17)
        assert result is not None
        assert isinstance(result, post.Mesh)
        assert len(result.node_ids) == 2
        assert len(result.element_ids) == 1


def test_meshes_select_solids_disjoint_from_shells_and_beams(meshes):
    """Test that solids, shells, and beams are disjoint selections."""
    solids = meshes.select_solids()
    shells = meshes.select_shells()
    beams = meshes.select_beams()

    solid_count = (
        0 if solids is None else (1 if isinstance(solids, post.Mesh) else len(solids))
    )
    shell_count = (
        0 if shells is None else (1 if isinstance(shells, post.Mesh) else len(shells))
    )
    beam_count = (
        0 if beams is None else (1 if isinstance(beams, post.Mesh) else len(beams))
    )

    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_12_0:
        # shells=5, solids=2, beams=1 → 8 out of 16 total
        assert solid_count == 2
        assert shell_count == 5
        assert beam_count == 1
    else:
        # shells=7, solids=2, beams=5 → 14 out of 18 total
        assert solid_count == 2
        assert shell_count == 7
        assert beam_count == 5

    # The sum of shape-specific meshes should not exceed total meshes
    assert solid_count + shell_count + beam_count <= len(meshes)


def test_meshes_select_solids_consistency_with_select(meshes):
    """Test that select_solids matches select with the solid elshape value."""
    from ansys.dpf.core import elements
    from ansys.dpf.core.check_version import get_server_version, meets_version

    from ansys.dpf import core

    if meets_version(get_server_version(core._global_server()), "12.0"):
        solid_value = elements._element_shapes.SOLID.value
    else:
        solid_value = elements._element_shapes_legacy.SOLID.value

    result_select = meshes.select(elshape=solid_value)
    result_select_solids = meshes.select_solids()

    if result_select is None:
        assert result_select_solids is None
    elif isinstance(result_select, post.Mesh):
        assert isinstance(result_select_solids, post.Mesh)
        assert len(result_select.node_ids) == len(result_select_solids.node_ids)
    else:
        assert isinstance(result_select_solids, post.Meshes)
        assert len(result_select) == len(result_select_solids)


def test_meshes_select_shells_consistency_with_select(meshes):
    """Test that select_shells matches select with the shell elshape value."""
    from ansys.dpf.core import elements
    from ansys.dpf.core.check_version import get_server_version, meets_version

    from ansys.dpf import core

    if meets_version(get_server_version(core._global_server()), "12.0"):
        shell_value = elements._element_shapes.SHELL.value
    else:
        shell_value = elements._element_shapes_legacy.SHELL.value

    result_select = meshes.select(elshape=shell_value)
    result_select_shells = meshes.select_shells()

    if result_select is None:
        assert result_select_shells is None
    elif isinstance(result_select, post.Mesh):
        assert isinstance(result_select_shells, post.Mesh)
        assert len(result_select.node_ids) == len(result_select_shells.node_ids)
    else:
        assert isinstance(result_select_shells, post.Meshes)
        assert len(result_select) == len(result_select_shells)


def test_meshes_select_beams_consistency_with_select(meshes):
    """Test that select_beams matches select with the beam elshape value."""
    from ansys.dpf.core import elements
    from ansys.dpf.core.check_version import get_server_version, meets_version

    from ansys.dpf import core

    if meets_version(get_server_version(core._global_server()), "12.0"):
        beam_value = elements._element_shapes.BEAM.value
    else:
        beam_value = elements._element_shapes_legacy.BEAM.value

    result_select = meshes.select(elshape=beam_value)
    result_select_beams = meshes.select_beams()

    if result_select is None:
        assert result_select_beams is None
    elif isinstance(result_select, post.Mesh):
        assert isinstance(result_select_beams, post.Mesh)
        assert len(result_select.node_ids) == len(result_select_beams.node_ids)
    else:
        assert isinstance(result_select_beams, post.Meshes)
        assert len(result_select) == len(result_select_beams)
