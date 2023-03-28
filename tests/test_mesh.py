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


def test_mesh_available_named_selections(mesh):
    ns = mesh.available_named_selections
    assert len(ns) == 1
    assert all([isinstance(n, str) for n in ns])


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
