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
            "mass_flow_rate",
            "static_pressure",
            "mean_static_pressure",
            "rms_static_pressure",
            "surface_heat_rate",
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
        result = getattr(fluent_simulation, result_name)(phases=[1])
        assert isinstance(result, post.DataFrame)
        result = getattr(fluent_simulation, result_name)(phases=["phase-1"])
        assert isinstance(result, post.DataFrame)
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
        result = getattr(cfx_simulation, result_name)(phases=[1])
        assert isinstance(result, post.DataFrame)
        result = getattr(cfx_simulation, result_name)(phases=[2])
        assert isinstance(result, post.DataFrame)

    def test_fluid_simulation_zones(self, fluent_simulation):
        from ansys.dpf.post.zone import Zones

        assert isinstance(fluent_simulation.zones, Zones)

    def test_fluid_simulation_species(self, fluent_simulation):
        from ansys.dpf.post.species import SpeciesDict

        assert isinstance(fluent_simulation.species, SpeciesDict)

    def test_fluid_simulation_phases(self, fluent_simulation):
        from ansys.dpf.post.phase import PhasesDict

        assert isinstance(fluent_simulation.phases, PhasesDict)

    def test_dataframe_plot_empty(self, fluent_simulation):
        result = fluent_simulation.wall_shear_stress()
        with pytest.raises(ValueError, match="No data to plot."):
            result.plot()

    def test_results_fluent_averaging_from_elemental(self, fluent_simulation):
        print(fluent_simulation)
        # ######## Elemental Result #################
        # Request on None
        result = fluent_simulation.enthalpy()
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental

        # Request on nodes
        result = fluent_simulation.enthalpy(location=post.locations.nodal)
        assert result.index.mesh_index.location == post.locations.nodal
        assert result._core_object[0].location == post.locations.nodal
        result = fluent_simulation.enthalpy_on_nodes()
        assert result.index.mesh_index.location == post.locations.nodal
        assert result._core_object[0].location == post.locations.nodal

        # Request on faces
        with pytest.raises(
            ValueError, match="Cannot query elemental results on faces."
        ):
            _ = fluent_simulation.enthalpy(location=post.locations.faces)

        # Request on cells
        result = fluent_simulation.enthalpy(location=post.locations.elemental)
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental
        result = fluent_simulation.enthalpy_on_cells()
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental

    def test_results_fluent_averaging_from_elemental_faces(self, fluent_simulation):
        print(fluent_simulation)
        # ######## ElementalFaces Result #################
        # Request on None
        result = fluent_simulation.static_pressure()
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental

        # Request on nodes
        result = fluent_simulation.static_pressure(location=post.locations.nodal)
        assert result.index.mesh_index.location == post.locations.nodal
        assert result._core_object[0].location == post.locations.nodal
        result = fluent_simulation.static_pressure_on_nodes()
        assert result.index.mesh_index.location == post.locations.nodal
        assert result._core_object[0].location == post.locations.nodal

        # Request on faces (requires filter-out of cell zones)
        # result = fluent_simulation.static_pressure(location=post.locations.faces)
        # assert result.index.mesh_index.location == post.locations.faces
        # assert result._core_object[0].location == post.locations.faces
        # result = fluent_simulation.static_pressure_on_faces()
        # assert result.index.mesh_index.location == post.locations.faces
        # assert result._core_object[0].location == post.locations.faces

        # Request on cells (requires filter-out of face zones)
        result = fluent_simulation.static_pressure(location=post.locations.elemental)
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental
        result = fluent_simulation.static_pressure_on_cells()
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental
