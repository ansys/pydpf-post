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

    @fixture
    def cfx_simulation(self):
        fluid_example_files = examples.download_cfx_heating_coil()
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

    def test_simulation_init(self, fluent_simulation):
        assert fluent_simulation is not None

    @pytest.mark.parametrize(
        "result_name",
        [
            "enthalpy",
            # "mass_flow_rate",
            "static_pressure",
            "mean_static_pressure",
            "rms_static_pressure",
            # "surface_heat_rate",
            "density",
            "temperature",
            "mean_temperature",
            "rms_temperature",
            "velocity",
            "mean_velocity",
            "rms_velocity",
        ],
    )
    def test_results_fluent(self, fluent_simulation, result_name):
        result = getattr(fluent_simulation, result_name)()
        assert isinstance(result, post.DataFrame)
        # result = getattr(fluent_simulation, result_name)(phases=[1])
        # assert isinstance(result, post.DataFrame)
        # result = getattr(fluent_simulation, result_name)(phases=["phase-1"])
        # assert isinstance(result, post.DataFrame)
        with pytest.raises(ValueError, match="is not a valid Phase ID or Phase name"):
            _ = getattr(fluent_simulation, result_name)(phases=[2])

    @pytest.mark.parametrize(
        "result_name",
        [
            "specific_heat",
            "epsilon",
            "enthalpy",
            "turbulent_kinetic_energy",
            "thermal_conductivity",
            "dynamic_viscosity",
            "turbulent_viscosity",
            "static_pressure",
            "total_pressure",
            "density",
            "entropy",
            "wall_shear_stress",
            "temperature",
            "total_temperature",
            "velocity",
        ],
    )
    def test_results_cfx(self, cfx_simulation, result_name):
        result = getattr(cfx_simulation, result_name)()
        assert isinstance(result, post.DataFrame)
        # result = getattr(cfx_simulation, result_name)(phases=[1])
        # assert isinstance(result, post.DataFrame)
        # result = getattr(cfx_simulation, result_name)(phases=[2])
        # assert isinstance(result, post.DataFrame)

    def test_fluid_simulation_zones(self, fluent_simulation):
        from ansys.dpf.post.zone import Zones

        assert isinstance(fluent_simulation.zones, Zones)

    def test_fluid_simulation_species(self, fluent_simulation):
        from ansys.dpf.post.species import SpeciesList

        assert isinstance(fluent_simulation.species, SpeciesList)

    def test_fluid_simulation_phases(self, fluent_simulation):
        from ansys.dpf.post.phase import Phases

        assert isinstance(fluent_simulation.phases, Phases)
