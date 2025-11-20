# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
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

import os

from ansys.dpf.core import examples
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_10_0,
)


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

    def test_simulation_str(self, fluent_simulation):
        assert fluent_simulation is not None
        assert str(fluent_simulation)

    @pytest.mark.skipif(
        (not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_10_0) and (os.name == "posix"),
        reason="Bug due to gatebin incompatibilities for servers <26.1",
    )
    def test_simulation_flprj(self):
        simulation = post.FluidSimulation(
            flprj=examples.download_fluent_axial_comp()["flprj"]
        )
        assert len(simulation.time_freq_support.time_frequencies) == 3
        result = simulation.static_pressure(all_sets=True)
        assert len(result.columns.set_ids) == 3

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

    def test_fluid_simulation_species(self, fluent_simulation):
        from ansys.dpf.post.species import SpeciesDict

        assert isinstance(fluent_simulation.species, SpeciesDict)

    def test_fluid_simulation_phases(self, fluent_simulation):
        from ansys.dpf.post.phase import PhasesDict

        assert isinstance(fluent_simulation.phases, PhasesDict)

    def test_fluid_simulation_result_unavailable(self, fluent_simulation):
        # print(fluent_simulation)
        with pytest.raises(ValueError, match="is not available."):
            _ = fluent_simulation.mass_fraction()

    def test_results_fluent_averaging_from_elemental(self, fluent_simulation):
        # print(fluent_simulation)
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
        # print(fluent_simulation)
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
        if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            with pytest.raises(
                ValueError,
                match="Querying an ElementalAndFaces result on "
                "faces currently requires the use of face zone ids",
            ):
                _ = fluent_simulation.static_pressure(location=post.locations.faces)
            with pytest.raises(
                ValueError,
                match="Querying an ElementalAndFaces result on "
                "faces currently requires the use of face zone ids",
            ):
                _ = fluent_simulation.static_pressure_on_faces()
        else:
            result = fluent_simulation.static_pressure(location=post.locations.faces)
            # print(result)
            assert result.index.mesh_index.location == post.locations.faces
            assert result._core_object[0].location == post.locations.faces
            # result._fc[0].plot()

            result = fluent_simulation.static_pressure_on_faces()
            # print(result)
            assert result.index.mesh_index.location == post.locations.faces
            assert result._core_object[0].location == post.locations.faces
            # result.plot()

        # Request on cells (requires filter-out of face zones)
        result = fluent_simulation.static_pressure(location=post.locations.elemental)
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental
        result = fluent_simulation.static_pressure_on_cells()
        assert result.index.mesh_index.location == "cells"
        assert result._core_object[0].location == post.locations.elemental

    def test_results_cfx_cross_locations_on_nodes(self, cfx_simulation):
        result = cfx_simulation.temperature_on_nodes(
            node_ids=cfx_simulation.mesh.node_ids
        )
        assert result.index.mesh_index.location == post.locations.nodal
        assert len(result.index.mesh_index) == cfx_simulation.mesh_info.num_nodes
        ref = """
  results          TEMP (K)           
  set_ids                 1           
    phase Water at 25 C (2) Copper (3)
 node_ids                             
        1        3.0533e+02           
        2        3.0389e+02           
        3        3.0155e+02           
        4        3.0201e+02           
        5        3.0340e+02           
        6        3.0487e+02           
      ...               ...        ...
"""  # noqa: W291, E501
        assert str(result) == ref
        result = cfx_simulation.temperature_on_nodes(
            cell_ids=cfx_simulation.mesh.element_ids
        )
        assert result.index.mesh_index.location == post.locations.nodal
        assert len(result.index.mesh_index) == cfx_simulation.mesh_info.num_nodes
        ref = """
  results          TEMP (K)           
  set_ids                 1           
    phase Water at 25 C (2) Copper (3)
 node_ids                             
     3149        3.0188e+02           
     4143        3.0039e+02           
     3140        3.0099e+02           
     3158        3.0127e+02           
     3154        3.0498e+02           
     4146        3.0049e+02           
      ...               ...        ...
"""  # noqa: W291, E501
        assert str(result) == ref
        result = cfx_simulation.temperature_on_nodes(
            face_ids=cfx_simulation.mesh.face_ids
        )
        assert result.index.mesh_index.location == post.locations.nodal
        ref = """
  results          TEMP (K)           
  set_ids                 1           
    phase Water at 25 C (2) Copper (3)
 node_ids                             
"""  # noqa: W291, E501
        assert str(result) == ref

    def test_results_cfx_cross_locations_on_faces(self, cfx_simulation):
        #         result = cfx_simulation.density_on_faces(
        #             cell_ids=cfx_simulation.mesh.element_ids
        #         )
        #         assert result.index.mesh_index.location == post.locations.faces
        #         assert result._fc[0].location == post.locations.faces
        #         ref = """
        #   results     RHO (kg*m^-3)
        #   set_ids                 1
        #     phase Water at 25 C (2) Copper (3)
        #  face_ids
        #         1        9.9700e+02
        #         2        9.9700e+02
        #         3        9.9700e+02
        #         4        9.9700e+02
        #         5        9.9700e+02
        #         6        9.9700e+02
        #       ...               ...        ...
        # """  # noqa: W291, E501
        #         assert str(result) == ref
        #         result.plot()
        result = cfx_simulation.temperature_on_faces(
            face_ids=cfx_simulation.mesh.face_ids
        )
        assert result.index.mesh_index.location == post.locations.faces
        ref = """
  results          TEMP (K)           
  set_ids                 1           
    phase Water at 25 C (2) Copper (3)
 face_ids                             
"""  # noqa: W291, E501
        assert str(result) == ref
        with pytest.raises(
            ValueError, match="Cannot plot a Dataframe with an empty mesh index."
        ):
            result.plot()

    def test_results_cfx_cross_locations_on_cells(self, cfx_simulation):
        result = cfx_simulation.temperature_on_cells(
            cell_ids=cfx_simulation.mesh.element_ids
        )
        assert result.index.mesh_index.location == "cells"
        assert len(result.index.mesh_index) == cfx_simulation.mesh_info.num_cells
        ref = """
  results          TEMP (K)           
  set_ids                 1           
    phase Water at 25 C (2) Copper (3)
 cell_ids                             
        1        3.0113e+02           
        2        3.0206e+02           
        3        3.0182e+02           
        4        3.0093e+02           
        5        3.0291e+02           
        6        3.0055e+02           
      ...               ...        ...
"""  # noqa: W291, E501
        assert str(result) == ref

    def test_results_fluent_cross_locations_on_nodes(self, fluent_simulation):
        result = fluent_simulation.density_on_nodes(
            node_ids=fluent_simulation.mesh.node_ids
        )
        assert result.index.mesh_index.location == post.locations.nodal
        ref = """
  results RHO (kg*m^-3)
  set_ids             1
 node_ids              
        1    1.0742e+00
        2    1.0436e+00
        3    1.0131e+00
        4    1.0327e+00
        5    1.0247e+00
        6    1.0445e+00
      ...           ...
"""  # noqa: W291, E501
        assert str(result) == ref
        result = fluent_simulation.density_on_nodes(
            cell_ids=fluent_simulation.mesh.element_ids
        )
        assert result.index.mesh_index.location == post.locations.nodal
        ref = """
  results RHO (kg*m^-3)
  set_ids             1
 node_ids              
      996    1.1041e+00
      894    1.1035e+00
      795    1.0982e+00
      903    1.0985e+00
      997    1.1097e+00
      895    1.1091e+00
      ...           ...
"""  # noqa: W291, E501
        assert str(result) == ref
        result = fluent_simulation.density_on_nodes(
            face_ids=fluent_simulation.mesh.face_ids
        )
        assert result.index.mesh_index.location == post.locations.nodal
        ref = """
  results RHO (kg*m^-3)
  set_ids             1
 node_ids              
    11291    1.3590e+00
    11416    1.3262e+00
    11455    1.3104e+00
    11325    1.3470e+00
    11348    1.2896e+00
    11388    1.2771e+00
      ...           ...
"""  # noqa: W291, E501
        assert str(result) == ref

    def test_results_fluent_cross_locations_on_faces(self, fluent_simulation):
        # TODO investigate wrong plot, wrong mesh index for dataframes
        # print(fluent_simulation)
        # with pytest.raises(
        #     ValueError,
        #     match="Querying an ElementalAndFaces result on "
        #     "faces currently requires the use of face zone ids",
        # ):
        #     _ = fluent_simulation.density_on_faces(
        #         cell_ids=fluent_simulation.mesh.element_ids
        #     )
        #         print(result)
        #         assert result.index.mesh_index.location == post.locations.faces
        #         # assert len(result.index.mesh_index.values) == fluent_simulation.mesh.num_faces
        #         ref = """
        #   results RHO (kg*m^-3)
        #   set_ids             1
        #  face_ids
        #         1    1.1095e+00
        #         2    1.1087e+00
        #         3    1.1098e+00
        #         4    1.0977e+00
        #         5    1.0949e+00
        #         6    1.1077e+00
        #       ...           ...
        # """  # noqa: W291, E501
        #         assert str(result) == ref
        #         result.plot()
        if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            with pytest.raises(
                ValueError,
                match="Querying an ElementalAndFaces result on "
                "faces currently requires the use of face zone ids",
            ):
                _ = fluent_simulation.density_on_faces(
                    face_ids=fluent_simulation.mesh.face_ids
                )
        else:
            result = fluent_simulation.density_on_faces(
                face_ids=fluent_simulation.mesh.face_ids
            )
            print(result)
            assert result.index.mesh_index.location == post.locations.faces
            # assert (
            #     len(result.index.mesh_index.values) == len(fluent_simulation.mesh.face_ids)
            # )  # TODO: why does this fail? Is the result not defined everywhere?
            ref = """
  results RHO (kg*m^-3)
  set_ids             1
 face_ids              
    39897    1.1077e+00
    39898    1.0892e+00
    39899    1.0821e+00
    39900    1.0761e+00
    39901    1.0721e+00
    39902    1.0724e+00
      ...           ...
"""  # noqa: W291, E501
            assert str(result) == ref

    #         result.plot()

    def test_results_fluent_cross_locations_on_cells(self, fluent_simulation):
        result = fluent_simulation.density_on_cells(
            cell_ids=fluent_simulation.mesh.element_ids
        )
        assert result.index.mesh_index.location == "cells"
        assert (
            len(result.index.mesh_index.values) == fluent_simulation.mesh.num_elements
        )
        ref = """
  results RHO (kg*m^-3)
  set_ids             1
 cell_ids              
        1    1.1095e+00
        2    1.1087e+00
        3    1.1098e+00
        4    1.0977e+00
        5    1.0949e+00
        6    1.1077e+00
      ...           ...
"""  # noqa: W291, E501
        assert str(result) == ref
        result.plot()

    def test_plot_result_on_zones(self, fluent_simulation):
        temperature = fluent_simulation.temperature(
            zone_ids=list(fluent_simulation.cell_zones.keys())
        )
        temperature.plot()

        temperature = fluent_simulation.temperature(
            qualifiers={"zone": list(fluent_simulation.cell_zones.keys())}
        )
        temperature.plot()

        temperature = fluent_simulation.temperature(
            zone_ids=list(fluent_simulation.face_zones.keys())
        )
        temperature.plot()

        temperature = fluent_simulation.temperature(
            qualifiers={"zone": list(fluent_simulation.face_zones.keys())}
        )
        temperature.plot()
