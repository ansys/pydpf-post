from ansys.dpf.core import examples
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
class TestFluidSimulation:
    @fixture
    def fluid_simulation(self):
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

    def test_simulation_init(self, fluid_simulation):
        assert fluid_simulation is not None

    @pytest.mark.parametrize(
        "result_name",
        [
            "enthalpy",
            "mass_flow_rate",
            "static_pressure",
            "mean_static_pressure",
            "rms_static_pressure",
            # "surface_heat_rate",  # Wait for fix
            "density",
            "temperature",
            "mean_temperature",
            "rms_temperature",
            "velocity",
            "mean_velocity",
            "rms_velocity",
        ],
    )
    def test_results(self, fluid_simulation, result_name):
        result = getattr(fluid_simulation, result_name)()
        assert result is not None
        assert isinstance(result, post.DataFrame)

    def test_fluid_simulation_zones(self, fluid_simulation):
        from ansys.dpf.post.zone import Zones

        assert isinstance(fluid_simulation.zones, Zones)

    def test_fluid_simulation_species(self, fluid_simulation):
        from ansys.dpf.post.species import SpeciesList

        assert isinstance(fluid_simulation.species, SpeciesList)

    def test_fluid_simulation_phases(self, fluid_simulation):
        from ansys.dpf.post.phase import Phases

        assert isinstance(fluid_simulation.phases, Phases)
