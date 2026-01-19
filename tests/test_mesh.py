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
import numpy as np
import pytest
from pytest import fixture

from ansys.dpf.post import FluidSimulation, Mesh, StaticMechanicalSimulation
from ansys.dpf.post.connectivity import ConnectivityListByIndex
from ansys.dpf.post.faces import Face
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_1,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_11_0,
)


@fixture
def mesh(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.mesh


@fixture
def fluent_mesh(fluid_fluent_elbow_steady_state):
    ds = dpf.DataSources()
    ds.set_result_file_path(
        fluid_fluent_elbow_steady_state["cas"][0],
        key="cas",
    )
    ds.add_file_path(
        fluid_fluent_elbow_steady_state["dat"][0],
        key="dat",
    )
    simulation = FluidSimulation(ds)
    return simulation.mesh  # noqa


def test_mesh_core_object(mesh):
    assert isinstance(mesh._core_object, dpf.MeshedRegion)
    assert mesh._core_object.nodes.n_nodes == 81
    with pytest.raises(ValueError, match="Tried to instantiate an empty Mesh."):
        _ = Mesh(None)


def test_mesh_node_ids(mesh):
    n_ids = mesh.node_ids
    assert len(n_ids) == 81
    assert all([isinstance(i, (int, np.integer)) for i in n_ids])


def test_mesh_element_ids(mesh):
    e_ids = mesh.element_ids
    assert len(e_ids) == 8
    assert all([isinstance(i, (int, np.integer)) for i in e_ids])


def test_mesh_num(mesh):
    assert mesh.num_nodes == 81
    assert mesh.num_elements == 8


def test_mesh_named_selections(mesh):
    ns = mesh.named_selections
    assert ns.keys()[0] == "_FIXEDSU"
    n = ns[ns.keys()[0]]
    assert len(n.ids) == 21


def test_mesh_unit(mesh):
    assert mesh.unit == "m"
    mesh.unit = "kg"
    assert mesh.unit == "kg"


def test_mesh_nodes(mesh):
    node_idx_0 = mesh.nodes[0]
    node_id_81 = mesh.nodes.by_id[81]

    assert len(mesh.nodes) == 81
    assert node_idx_0.index == 0
    assert node_id_81.id == 81
    assert all(np.isclose(node_idx_0.coordinates, [0.015, 0.045, 0.015]))
    assert all(np.isclose(node_id_81.coordinates, [0.03, 0.045, 0.0075]))
    assert np.array_equal([node.index for node in mesh.nodes], list(range(81)))
    assert mesh.get_node_by_id(node_idx_0.id).index == 0


def test_mesh_elements(mesh):
    elem_idx_0 = mesh.elements[0]
    elem_id_8 = mesh.elements.by_id[8]

    assert len(mesh.elements) == 8
    assert elem_idx_0.index == 0
    assert elem_id_8.id == 8
    assert elem_idx_0.num_nodes == 20
    assert elem_idx_0.type == dpf.element_types.Hex20
    assert elem_idx_0.shape == "solid"
    assert mesh.get_element_by_id(elem_idx_0.id).index == 0


def test_mesh_nodal_connectivity(mesh):
    conn_node_get_idx = mesh.node_to_element_connectivity
    conn_node_get_id = mesh.node_to_element_ids_connectivity

    node_idx_0 = mesh.nodes[0]
    conn_idx = node_idx_0.to_element_connectivity
    conn_id = [mesh.elements[idx].id for idx in conn_idx]

    assert len(conn_idx) == len(conn_id)
    assert np.array_equal(conn_idx, conn_node_get_idx[0])
    assert np.array_equal(conn_id, conn_node_get_id[0])


def test_mesh_elemental_connectivity(mesh):
    conn_elem_get_idx = mesh.element_to_node_connectivity
    conn_elem_get_id = mesh.element_to_node_ids_connectivity

    elem_idx_0 = mesh.elements[0]
    conn_idx = elem_idx_0.to_node_connectivity
    conn_id = [mesh.nodes[idx].id for idx in conn_idx]

    assert len(conn_idx) == len(conn_id)
    assert np.array_equal(conn_idx, conn_elem_get_idx[0])
    assert np.array_equal(conn_id, conn_elem_get_id[0])


def test_mesh_str(mesh):
    txt = str(mesh)
    assert (
        txt
        == "DPF  Mesh: \n  81 nodes \n  8 elements \n  Unit: m \n  With solid (3D) elements"
    )


def test_mesh_coordinates(mesh):
    coord = mesh.coordinates
    # print(coord)
    ref = """
             results  coord (m)
 node_ids components           
        1          X 1.5000e-02
                   Y 4.5000e-02
                   Z 1.5000e-02
        2          X 1.5000e-02
                   Y 4.5000e-02
                   Z 0.0000e+00
      ...        ...        ...
"""  # noqa
    assert str(coord) == ref


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_1,
    reason="PropertyFieldsContainer supported for server versions greater than 8.1",
)
def test_mesh_materials(mesh):
    materials = mesh.materials
    ref = """
     results material_id
 element_ids            
           5           1
           6           1
           1           1
           2           1
           7           1
           8           1
         ...         ...
"""  # noqa
    assert str(materials) == ref
    materials_5 = materials.select(element_ids=[5])
    ref = """
     results material_id
 element_ids            
           5           1
"""  # noqa
    assert str(materials_5) == ref


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_1,
    reason="PropertyFieldsContainer supported for server versions greater than 8.1",
)
def test_mesh_element_types(mesh):
    element_types = mesh.element_types
    # print(element_types)
    ref = """
     results elem_type_id
 element_ids             
           5            1
           6            1
           1            1
           2            1
           7            1
           8            1
         ...          ...
"""  # noqa
    assert str(element_types) == ref


def test_mesh_plot(mesh):
    mesh.plot()


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Available starting DPF 7.0",
)
def test_mesh_faces(fluent_mesh):
    assert fluent_mesh.num_faces == 94288
    assert len(fluent_mesh.face_ids) == 94288
    assert len(fluent_mesh.faces) == 94288
    assert isinstance(fluent_mesh.faces[0], Face)
    first_face = fluent_mesh.get_face_by_id(fluent_mesh.face_ids[0])
    assert isinstance(first_face, Face)
    assert first_face.index == 0
    assert first_face.id == fluent_mesh.face_ids[0]
    assert isinstance(fluent_mesh.face_to_node_connectivity, ConnectivityListByIndex)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_11_0:
        assert fluent_mesh.face_to_node_connectivity[0] == [2921, 25, 20]
        assert fluent_mesh.face_to_node_ids_connectivity[1] == [23, 2922, 21]
    else:
        assert fluent_mesh.face_to_node_connectivity[0] == [20, 25, 2921]
        assert fluent_mesh.face_to_node_ids_connectivity[1] == [21, 2922, 23]
    assert isinstance(
        fluent_mesh.face_to_node_ids_connectivity, ConnectivityListByIndex
    )
