from ansys.dpf.core import examples
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_0
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_0,
    # TODO: change to SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
class TestFluidSimulation:
    @fixture
    def fluid_simulation(self):
        fluid_example_files = examples.download_fluent_files()
        ds = dpf.DataSources()
        ds.set_result_file_path(
            fluid_example_files["cas"],
            key="cas",
        )
        ds.add_file_path(
            fluid_example_files["dat"],
            key="dat",
        )
        return post.FluidSimulation(ds)  # noqa

    def test_simulation_init(self, fluid_simulation):
        assert fluid_simulation is not None

    def test_density(self, fluid_simulation):
        result = fluid_simulation.density()
        assert result is not None
        print(result)
