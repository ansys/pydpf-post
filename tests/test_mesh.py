import ansys.dpf.core as core
from pytest import fixture

from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation


@fixture
def mesh(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.mesh


def test_mesh_core_object(mesh):
    assert isinstance(mesh._core_object, core.MeshedRegion)
    assert mesh._core_object.nodes.n_nodes == 81
