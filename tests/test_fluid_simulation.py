from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post.common import AvailableSimulationTypes


@fixture
def static_simulation(static_rst):
    return post.load_simulation(
        data_sources=static_rst,
        simulation_type=AvailableSimulationTypes.steady_fluid,
    )


def test_simulation_init(static_rst):
    simulation = post.FluidSimulation(static_rst)
    assert simulation is not None
