import ansys.dpf.core as dpf
import numpy as np
from pytest import fixture

from ansys.dpf.post import StaticMechanicalSimulation


@fixture
def mesh(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.mesh


def test_mesh_core_object(mesh):
    assert isinstance(mesh._core_object, dpf.MeshedRegion)
    assert mesh._core_object.nodes.n_nodes == 81


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
    ns = mesh.named_selections.keys()
    assert len(ns) == 1
    assert all([isinstance(n, str) for n in ns])
    assert ns[0] == "_FIXEDSU"
    assert len(mesh.named_selections[ns[0]].ids) == 21


def test_mesh_unit(mesh):
    assert mesh.unit == "m"


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
    assert elem_idx_0.type == dpf.element_types.Hex20.value
    assert elem_idx_0.shape == "solid"
    assert mesh.get_element_by_id(elem_idx_0.id).index == 0


def test_mesh_nodal_connectivity(mesh):
    conn_node_get_idx = mesh.node_to_element_connectivity
    conn_node_get_id = mesh.node_to_element_ids_connectivity

    node_idx_0 = mesh.nodes[0]
    conn_idx = node_idx_0.nodal_connectivity
    conn_id = [mesh.elements[idx].id for idx in conn_idx]

    assert len(conn_idx) == len(conn_id)
    assert np.array_equal(conn_idx, conn_node_get_idx[0])
    assert np.array_equal(conn_id, conn_node_get_id[0])


def test_mesh_elemental_connectivity(mesh):
    conn_elem_get_idx = mesh.element_to_node_connectivity
    conn_elem_get_id = mesh.element_to_node_ids_connectivity

    elem_idx_0 = mesh.elements[0]
    conn_idx = elem_idx_0.connectivity
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
    print(coord)
    ref = """
                 results   coord (m)
    node_ids  components            
           1           X  1.5000e-02
                       Y  4.5000e-02
                       Z  1.5000e-02
           2           X  1.5000e-02
                       Y  4.5000e-02
                       Z  0.0000e+00
         ...
"""  # noqa
    assert str(coord) == ref
