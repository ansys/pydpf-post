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

import csv
import dataclasses
from itertools import groupby
import os.path
import pathlib
from typing import Optional, Union

import ansys.dpf.core as dpf
from ansys.dpf.core import (
    Field,
    FieldsContainer,
    MeshedRegion,
    Scoping,
    element_types,
    natures,
    operators,
    shell_layers,
)
from ansys.dpf.gate.common import locations
import numpy as np
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post import StaticMechanicalSimulation
from ansys.dpf.post.common import AvailableSimulationTypes, elemental_properties
from ansys.dpf.post.index import ref_labels
from ansys.dpf.post.meshes import Meshes
from ansys.dpf.post.result_workflows._component_helper import ResultCategory
from ansys.dpf.post.result_workflows._utils import (
    AveragingConfig,
    _CreateOperatorCallable,
)
from ansys.dpf.post.selection import _WfNames
from ansys.dpf.post.simulation import MechanicalSimulation, Simulation
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_2,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_10_0,
    ReferenceCsvFilesNodal,
)


def is_principal(mode: Optional[str]) -> bool:
    return mode == "principal"


def is_equivalent(mode: Optional[str]) -> bool:
    return mode == "equivalent"


def mode_suffix(mode: Optional[str]) -> str:
    if mode == "equivalent":
        return "_eqv_von_mises"
    elif mode == "principal":
        return "_principal"
    return ""


def get_expected_elemental_average_skin_value(
    element_id: int,
    solid_mesh: MeshedRegion,
    skin_mesh: MeshedRegion,
    elemental_nodal_solid_data_field: Field,
    is_principal_strain_result: bool,
) -> dict[int, float]:
    """
    Get the average skin value of all the skin elements that belong to the solid
    element with the element_id.
    Returns a dictionary with the expected average skin values indexed
    by the skin element id.
    """

    elemental_nodal_solid_data = elemental_nodal_solid_data_field.get_entity_data_by_id(
        element_id
    )

    # Find all nodes connected to this solid element
    connected_node_indices = (
        solid_mesh.elements.connectivities_field.get_entity_data_by_id(element_id)
    )
    connected_node_ids = solid_mesh.nodes.scoping.ids[connected_node_indices]

    # Find the skin elements attached to all the nodes of the solid element
    # Note: This includes elements that are not skin elements of
    # a particular solid element (because not all nodes are part of the solid element)
    skin_element_ids = set()
    for connected_node_id in connected_node_ids:
        skin_element_index = (
            skin_mesh.nodes.nodal_connectivity_field.get_entity_data_by_id(
                connected_node_id
            )
        )
        skin_element_ids.update(skin_mesh.elements.scoping.ids[skin_element_index])

    expected_average_skin_values = {}
    for skin_element_id in skin_element_ids:
        # Go through all the skin element candidates and check if all their nodes are
        # part of the solid element.
        skin_node_indices = (
            skin_mesh.elements.connectivities_field.get_entity_data_by_id(
                skin_element_id
            )
        )
        node_ids = skin_mesh.nodes.scoping.ids[skin_node_indices]

        node_missing = False
        for node_id in node_ids:
            if node_id not in connected_node_ids:
                node_missing = True
                break

        if node_missing:
            # If a node is missing the skin element does not belong to the solid element
            continue

        # Get the elementary indices of the connected nodes in the solid
        # ElementalNodal data
        indices_to_average = []
        for idx, node_id in enumerate(connected_node_ids):
            if node_id in node_ids:
                indices_to_average.append(idx)

        # Skip indices that are out of bounds (these are the mid side nodes)
        indices_to_average = [
            idx for idx in indices_to_average if idx < len(elemental_nodal_solid_data)
        ]
        if elemental_nodal_solid_data_field.component_count > 1:
            data_to_average = elemental_nodal_solid_data[indices_to_average, :]
        else:
            data_to_average = elemental_nodal_solid_data[indices_to_average]

        average = np.mean(data_to_average, axis=0)
        if is_principal_strain_result:
            # Workaround: The principal operator divides
            # offdiagonal components by 2 if the input field has
            # a integer "strain" property set. Since int
            # field properties are not exposed in python, division is done
            # here before the the data is passed to the principle operator
            average[3:7] = average[3:7] / 2

        expected_average_skin_values[skin_element_id] = average
    return expected_average_skin_values


def get_expected_skin_results(
    create_operator_callable: _CreateOperatorCallable,
    element_ids: list[int],
    elemental_nodal_fc: FieldsContainer,
    meshed_region: MeshedRegion,
    mode: Optional[str],
    result_name: str,
):
    expected_skin_values = {}
    for element_id in element_ids:
        expected_skin_values_per_element = get_expected_elemental_average_skin_value(
            element_id=element_id,
            solid_mesh=elemental_nodal_fc[0].meshed_region,
            skin_mesh=meshed_region,
            elemental_nodal_solid_data_field=elemental_nodal_fc[0],
            is_principal_strain_result=is_principal(mode)
            and result_name == "elastic_strain",
        )

        if is_principal(mode) or (
            is_equivalent(mode) and result_name != "elastic_strain"
        ):
            # We need to put the expected skin values in a Field, to compute
            # the equivalent or principal values with a dpf operator.
            # For the elastic strain result, the invariants are computed before the
            # averaging
            field = Field(nature=natures.symmatrix)
            for (
                skin_element_id,
                expected_skin_value,
            ) in expected_skin_values_per_element.items():
                field.append(list(expected_skin_value), skin_element_id)
            if is_principal(mode):
                invariant_op = operators.invariant.principal_invariants()
                invariant_op.inputs.field(field)
                field_out = invariant_op.outputs.field_eig_1()
            else:
                invariant_op = create_operator_callable(name="eqv")
                invariant_op.inputs.field(field)
                field_out = invariant_op.outputs.field()

            for skin_element_id in expected_skin_values_per_element:
                expected_skin_values_per_element[
                    skin_element_id
                ] = field_out.get_entity_data_by_id(skin_element_id)
        expected_skin_values.update(expected_skin_values_per_element)
    return expected_skin_values


def get_and_check_elemental_skin_results(
    static_simulation: StaticMechanicalSimulation,
    fc_elemental_nodal: FieldsContainer,
    result_name: str,
    mode: str,
    element_ids: list[int],
    skin: Union[list[int], bool],
    expand_cyclic: bool,
):
    """
    Get the elemental skin results and check if they match the
    expected average skin values.
    """

    # Not all the simulation types have the expand cyclic option
    kwargs = {}
    if expand_cyclic:
        kwargs["expand_cyclic"] = True
    result_skin_scoped_elemental = getattr(
        static_simulation, f"{result_name}{mode_suffix(mode)}_elemental"
    )(set_ids=[1], skin=skin, **kwargs)

    if is_equivalent(mode) and result_name == "elastic_strain":
        # For the elastic strain result, the equivalent strains are computed
        # before the averaging
        invariant_op = static_simulation._model.operator(name="eqv_fc")
        invariant_op.inputs.fields_container(fc_elemental_nodal)
        fc_elemental_nodal = invariant_op.outputs.fields_container()

    expected_skin_values = get_expected_skin_results(
        create_operator_callable=static_simulation._model.operator,
        element_ids=element_ids,
        elemental_nodal_fc=fc_elemental_nodal,
        meshed_region=result_skin_scoped_elemental._fc[0].meshed_region,
        mode=mode,
        result_name=result_name,
    )

    for skin_element_id, expected_skin_value in expected_skin_values.items():
        actual_skin_value = result_skin_scoped_elemental._fc[0].get_entity_data_by_id(
            skin_element_id
        )
        assert np.allclose(actual_skin_value, expected_skin_value)


operator_map = {"stress": "S", "elastic_strain": "EPEL", "displacement": "U"}


def get_elemental_nodal_results(
    simulation: Simulation,
    result_name: str,
    scoping: Scoping,
    mode: str,
    expand_cyclic: bool,
):
    time_id = 1
    if expand_cyclic:
        result_result_scoped_elemental = getattr(
            simulation, f"{result_name}{mode_suffix(mode)}"
        )(set_ids=[time_id], expand_cyclic=True)

        return result_result_scoped_elemental._fc

    else:
        elemental_nodal_result_op = simulation._model.operator(
            name=operator_map[result_name]
        )
        if scoping is not None:
            elemental_nodal_result_op.inputs.mesh_scoping(scoping)
        elemental_nodal_result_op.inputs.time_scoping([time_id])
        elemental_nodal_result_op.inputs.requested_location("ElementalNodal")

        return elemental_nodal_result_op.outputs.fields_container()


def get_expected_nodal_averaged_skin_results(
    skin_mesh: MeshedRegion,
    solid_elemental_nodal_results: Field,
    is_principal_strain_result: bool,
):
    nodal_averaged_skin_values = Field(
        nature=solid_elemental_nodal_results.field_definition.dimensionality.nature,
        location=locations.nodal,
    )
    solid_mesh = solid_elemental_nodal_results.meshed_region

    all_midside_node_ids = get_all_midside_node_ids(
        solid_elemental_nodal_results.meshed_region
    )
    for skin_node_id in skin_mesh.nodes.scoping.ids:
        if skin_node_id in all_midside_node_ids:
            # Skip midside nodes. We don't extract results for
            # midside nodes
            continue
        solid_elements_indices = (
            solid_mesh.nodes.nodal_connectivity_field.get_entity_data_by_id(
                skin_node_id
            )
        )
        solid_elements_ids = solid_mesh.elements.scoping.ids[solid_elements_indices]

        solid_elemental_nodal_value = {}

        # Get the correct elemental nodal value for each adjacent solid element
        # will be later used to average the nodal values
        for solid_element_id in solid_elements_ids:
            if solid_element_id not in solid_elemental_nodal_results.scoping.ids:
                # Solid element is connected to the node but does not have
                # a value because it was not selected with the scoping
                continue
            solid_element_data = solid_elemental_nodal_results.get_entity_data_by_id(
                solid_element_id
            )
            connected_node_indices = (
                solid_mesh.elements.connectivities_field.get_entity_data_by_id(
                    solid_element_id
                )
            )
            connected_node_ids = solid_mesh.nodes.scoping.ids[connected_node_indices]
            node_index_in_elemental_data = np.where(connected_node_ids == skin_node_id)[
                0
            ][0]
            solid_elemental_nodal_value[solid_element_id] = solid_element_data[
                node_index_in_elemental_data
            ]

        # Average over the adjacent skin elements
        values_to_average = np.empty(
            shape=(0, solid_elemental_nodal_results.component_count)
        )
        skin_element_indices = (
            skin_mesh.nodes.nodal_connectivity_field.get_entity_data_by_id(skin_node_id)
        )
        skin_element_ids = skin_mesh.elements.scoping.ids[skin_element_indices]
        for skin_element_id in skin_element_ids:
            matching_solid_element_id = None
            connected_skin_node_indices = (
                skin_mesh.elements.connectivities_field.get_entity_data_by_id(
                    skin_element_id
                )
            )
            connected_skin_node_ids = skin_mesh.nodes.scoping.ids[
                connected_skin_node_indices
            ]

            for solid_element_id in solid_elemental_nodal_value.keys():
                connected_solid_node_indices = (
                    solid_mesh.elements.connectivities_field.get_entity_data_by_id(
                        solid_element_id
                    )
                )
                connected_solid_node_ids = solid_mesh.nodes.scoping.ids[
                    connected_solid_node_indices
                ]

                if set(connected_skin_node_ids).issubset(connected_solid_node_ids):
                    matching_solid_element_id = solid_element_id
                    break

            if matching_solid_element_id is None:
                raise RuntimeError(
                    f"No matching solid element found for skin element {skin_element_id}"
                )
            skin_nodal_value = solid_elemental_nodal_value[matching_solid_element_id]
            values_to_average = np.vstack((values_to_average, skin_nodal_value))
        skin_values = np.mean(values_to_average, axis=0)
        if is_principal_strain_result:
            # Workaround: The principal operator divides
            # offdiagonal components by 2 if the input field has
            # a integer "strain" property set. Since int
            # field properties are not exposed in python, division is done
            # here before the the data is passed to the principle operator
            skin_values[3:7] = skin_values[3:7] / 2

        nodal_averaged_skin_values.append(list(skin_values), skin_node_id)
    return nodal_averaged_skin_values


def get_all_midside_node_ids(mesh: MeshedRegion):
    all_midside_nodes = set()
    for element in mesh.elements:
        element_descriptor = element_types.descriptor(element.type)

        all_node_ids = element.node_ids
        for idx, node_id in enumerate(all_node_ids):
            if idx >= element_descriptor.n_corner_nodes:
                all_midside_nodes.add(node_id)

    return all_midside_nodes


def get_expected_nodal_skin_results(
    simulation: MechanicalSimulation,
    result_name: str,
    mode: Optional[str],
    skin_mesh: MeshedRegion,
    expand_cyclic: bool,
    elemental_nodal_results: Optional[FieldsContainer] = None,
):
    # We have two options to get nodal data:
    # 1)    Request nodal location directly from the operator.
    #       This way we get the nodal data of the full mesh scoped to
    #       the elements in the element scope.
    # 2)    Get the elemental nodal data and then
    #       average it to nodal. We get different results at the boundaries
    #       of the element scope compared to 1). This is because the averaging cannot take into
    #       account the elemental nodal data outside of the element scope. Therefore, the
    #       averaged node data at the boundaries is different.
    # Currently, the skin workflow requests elemental nodal data and then averages it to nodal,
    # which corresponds to the case 2 above.

    # If we don't get a elemental_nodal_result_op, this means the result does not support
    # elemental nodal evaluation. Currently used only for displacement.
    # In this case we just get the nodal results directly (case 1 above)
    time_id = 1
    if elemental_nodal_results is None:
        assert result_name == "displacement"
        kwargs = {}
        if expand_cyclic:
            kwargs["expand_cyclic"] = True
        nodal_field = simulation.displacement(set_ids=[time_id], **kwargs)._fc[0]
    else:
        if is_equivalent(mode) and result_name == "elastic_strain":
            # For elastic strain results, the computation of the equivalent
            # value happens before the averaging.
            invariant_op = simulation._model.operator(name="eqv_fc")
            invariant_op.inputs.fields_container(elemental_nodal_results)
            fields_container = invariant_op.outputs.fields_container()
        else:
            fields_container = elemental_nodal_results

        nodal_field = get_expected_nodal_averaged_skin_results(
            skin_mesh=skin_mesh,
            solid_elemental_nodal_results=fields_container[0],
            is_principal_strain_result=is_principal(mode)
            and result_name == "elastic_strain",
        )

    if is_principal(mode):
        invariant_op = simulation._model.operator(name="invariants")
        invariant_op.inputs.field(nodal_field)
        nodal_field = invariant_op.outputs.field_eig_1()

    if is_equivalent(mode) and result_name != "elastic_strain":
        invariant_op = simulation._model.operator(name="eqv")
        invariant_op.inputs.field(nodal_field)
        nodal_field = invariant_op.outputs.field()
    return nodal_field


@fixture
def static_simulation(static_rst):
    return post.load_simulation(
        data_sources=static_rst,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )


@fixture
def mixed_shell_solid_simulation(mixed_shell_solid_model):
    return post.load_simulation(
        data_sources=mixed_shell_solid_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )


@fixture
def mixed_shell_solid_with_contact_simulation(mixed_shell_solid_with_contact_model):
    return post.load_simulation(
        data_sources=mixed_shell_solid_with_contact_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )


@fixture
def two_cubes_contact_simulation(two_cubes_contact_model):
    return post.load_simulation(
        data_sources=two_cubes_contact_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )


@fixture
def transient_simulation(plate_msup):
    return post.load_simulation(
        data_sources=plate_msup,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )


@fixture
def modal_simulation(modalallkindofcomplexity):
    return post.load_simulation(
        data_sources=modalallkindofcomplexity,
        simulation_type=AvailableSimulationTypes.modal_mechanical,
    )


@fixture
def harmonic_simulation(complex_model):
    return post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    )


@fixture
def cyclic_static_simulation(simple_cyclic):
    return post.StaticMechanicalSimulation(simple_cyclic)


def test_simulation_init(static_rst):
    simulation = post.StaticMechanicalSimulation(static_rst)
    assert simulation is not None
    simulation = post.TransientMechanicalSimulation(static_rst)
    assert simulation is not None
    simulation = post.ModalMechanicalSimulation(static_rst)
    assert simulation is not None
    simulation = post.HarmonicMechanicalSimulation(static_rst)
    assert simulation is not None


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
    reason="Available starting DPF 4.0",
)
def test_simulation_init_with_server(static_rst, grpc_server):
    simulation = post.StaticMechanicalSimulation(static_rst, server=grpc_server)
    assert simulation is not None
    assert simulation._model._server != dpf.SERVER
    assert simulation._model._server == grpc_server


def test_simulation_units(static_simulation):
    assert static_simulation._units is None
    assert static_simulation.units is not None
    assert static_simulation.units["time"] == "s"
    assert static_simulation.units["length"] == "m"


def test_simulation_results(static_simulation):
    results = static_simulation.results
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
        assert len(results) == 12
    else:
        assert len(results) == 13
    assert all(
        isinstance(x, dpf.result_info.available_result.AvailableResult) for x in results
    )


def test_simulation_geometries(static_simulation):
    geometries = static_simulation.geometries
    assert geometries == []


# def test_simulation_boundary_conditions(static_simulation):
#     boundary_conditions = static_simulation.boundary_conditions
#     assert boundary_conditions == []
#
#
# def test_simulation_loads(static_simulation):
#     loads = static_simulation.loads
#     assert loads == []


def test_simulation_mesh(static_simulation):
    mesh = static_simulation.mesh
    assert isinstance(mesh, post.mesh.Mesh)


def test_simulation_named_selections(static_simulation):
    named_selections = static_simulation.named_selections
    assert len(named_selections) == 1
    assert all(isinstance(x, str) for x in named_selections)


def test_simulation_active_selection(static_simulation):
    assert static_simulation.active_selection is None
    selection = post.selection.Selection()
    static_simulation.active_selection = selection
    assert static_simulation.active_selection == selection
    static_simulation.deactivate_selection()
    assert static_simulation.active_selection is None


def test_simulation_plot(static_simulation):
    static_simulation.plot(cpos="xy")


def test_simulation_split_mesh_by_properties(allkindofcomplexity):
    simulation = post.StaticMechanicalSimulation(allkindofcomplexity)
    meshes = simulation.split_mesh_by_properties(
        properties=[
            elemental_properties.material,
            elemental_properties.element_shape,
        ]
    )
    assert isinstance(meshes, Meshes)
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0:
        assert len(meshes) == 18
    else:
        assert len(meshes) == 16
    meshes = simulation.split_mesh_by_properties(
        properties={
            elemental_properties.material: 1,
            elemental_properties.element_shape: [0, 1],
        }
    )
    assert isinstance(meshes, Meshes)
    assert len(meshes) == 2
    meshes = simulation.split_mesh_by_properties(
        properties={
            elemental_properties.material: 1,
            elemental_properties.element_shape: [0, 2],
        }
    )
    assert isinstance(meshes, post.Mesh)
    meshes = simulation.split_mesh_by_properties(
        properties={
            elemental_properties.material: 22,
            elemental_properties.element_shape: [0, 2],
        }
    )
    assert meshes is None


class TestStaticMechanicalSimulation:
    def test_cyclic(self, simple_cyclic):
        simulation = post.StaticMechanicalSimulation(simple_cyclic)
        result = simulation.stress(expand_cyclic=False)
        # print(result)
        assert "base_sector" in result.columns.names
        result = simulation.stress(expand_cyclic=True)
        # print(result)
        assert "base_sector" not in result.columns.names

    def test_multi_stage(self, multi_stage_cyclic):
        simulation = post.StaticMechanicalSimulation(multi_stage_cyclic)
        result = simulation.stress(expand_cyclic=False)
        # print(result)
        assert "base_sector" in result.columns.names
        assert "stage" in result.columns.names
        result = simulation.stress(expand_cyclic=True)
        # print(result)
        assert "base_sector" not in result.columns.names
        assert "stage" not in result.columns.names

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
        reason="Available starting DPF 4.0",
    )
    def test_with_grpc_server(self, static_rst, grpc_server):
        simulation = post.StaticMechanicalSimulation(static_rst, server=grpc_server)
        assert simulation._model._server != dpf.SERVER
        _ = simulation.displacement()
        _ = simulation.displacement(skin=True)

    def test_times_argument(self, static_simulation):
        _ = static_simulation.displacement(times=1)
        _ = static_simulation.displacement(times=1.0)
        _ = static_simulation.displacement(times=[1])
        _ = static_simulation.displacement(times=[1.0])
        with pytest.raises(
            ValueError, match="Argument times must contain numeric values only."
        ):
            _ = static_simulation.displacement(times=[0.0, 1, "test"])
        with pytest.raises(
            TypeError, match="Argument times must be a number or a list of numbers."
        ):
            _ = static_simulation.displacement(times="test")

    def test_warning_empty(self, static_simulation):
        with pytest.warns(expected_warning=UserWarning, match="empty"):
            _ = static_simulation.displacement(
                components=1, node_ids=[1001, 1002, 1003]
            )

    def test_raise_mutually_exclusive(self, static_simulation):
        with pytest.raises(ValueError, match="exclusive"):
            _ = static_simulation.displacement(node_ids=[42], element_ids=[1])
        with pytest.raises(ValueError, match="exclusive"):
            _ = static_simulation.displacement(load_steps=[1], set_ids=[1])

    def test_displacement(self, static_simulation):
        displacement_x = static_simulation.displacement(
            components=["X"], node_ids=[42, 43, 44]
        )
        assert len(displacement_x._fc) == 1
        assert displacement_x._fc.get_time_scoping().ids == [1]
        field = displacement_x._fc[0]
        op = static_simulation._model.operator("UX")
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [42, 43, 44], server=static_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

        with pytest.raises(
            ValueError, match="Sub-step 2 of load-step 1 does not exist."
        ):
            _ = static_simulation.displacement(
                components=["2"],
                named_selections=static_simulation.named_selections[0],
                load_steps=(1, 2),
            )

        displacement_y = static_simulation.displacement(
            components=["2"],
            named_selections=static_simulation.named_selections[0],
            load_steps=(1, 1),
        )
        assert len(displacement_y._fc) == 1
        assert displacement_y._fc.get_time_scoping().ids == [1]
        field = displacement_y._fc[0]
        op = static_simulation._model.operator("UY")
        mesh_scoping = dpf.mesh_scoping_factory.named_selection_scoping(
            static_simulation.named_selections[0],
            server=static_simulation._model._server,
            model=static_simulation._model,
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (21,)
        assert np.allclose(field.data, field_ref.data)

        displacement_z = static_simulation.displacement(
            components="Z",
            named_selections=static_simulation.named_selections[0],
            load_steps=(1, 1),
        )
        assert len(displacement_z._fc) == 1
        assert displacement_z._fc.get_time_scoping().ids == [1]
        field = displacement_z._fc[0]
        op = static_simulation._model.operator("UZ")
        mesh_scoping = dpf.mesh_scoping_factory.named_selection_scoping(
            static_simulation.named_selections[0],
            server=static_simulation._model._server,
            model=static_simulation._model,
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (21,)
        assert np.allclose(field.data, field_ref.data)

        displacement_z = static_simulation.displacement(
            components="Z",
            element_ids=[1, 2, 3],
            set_ids=1,
        )
        assert len(displacement_z._fc) == 1
        assert displacement_z._fc.get_time_scoping().ids == [1]
        field = displacement_z._fc[0]
        op = static_simulation._model.operator("UZ")
        mesh_scoping = dpf.mesh_scoping_factory.elemental_scoping(
            element_ids=[1, 2, 3],
            server=static_simulation._model._server,
        )
        mesh_scoping = dpf.operators.scoping.transpose(
            mesh_scoping=mesh_scoping,
            meshed_region=static_simulation.mesh._meshed_region,
            inclusive=1,
        ).eval()
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field_ref.data.shape == (44,)
        assert field.component_count == 1
        assert field.data.shape == (44,)
        assert np.allclose(field.data, field_ref.data)

        displacement_norm = static_simulation.displacement(
            norm=True, node_ids=[42, 43, 44], set_ids=[1]
        )
        assert len(displacement_norm._fc) == 1
        assert displacement_norm._fc.get_time_scoping().ids == [1]
        field = displacement_norm._fc[0]
        op = static_simulation._model.operator("U")
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [42, 43, 44], server=static_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        norm_op = static_simulation._model.operator("norm_fc")
        norm_op.connect(0, op.outputs.fields_container)
        field_ref = norm_op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress(self, static_simulation):
        stress_x = static_simulation.stress(components=1)
        assert len(stress_x._fc) == 1
        assert stress_x._fc.get_time_scoping().ids == [1]
        field = stress_x._fc[0]
        op = static_simulation._model.operator("SX")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (64,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_elemental(self, static_simulation):
        stress_x = static_simulation.stress_elemental(components=1)
        assert len(stress_x._fc) == 1
        assert stress_x._fc.get_time_scoping().ids == [1]
        field = stress_x._fc[0]
        op = static_simulation._model.operator("SX")
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (8,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_nodal(self, static_simulation):
        stress_x = static_simulation.stress_nodal(components=1)
        assert len(stress_x._fc) == 1
        assert stress_x._fc.get_time_scoping().ids == [1]
        field = stress_x._fc[0]
        op = static_simulation._model.operator("SX")
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (81,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal(self, static_simulation):
        result = static_simulation.stress_principal(components=1)
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("S1")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (64,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_nodal(self, static_simulation):
        result = static_simulation.stress_principal_nodal(components=2)
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("S2")
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (81,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_elemental(self, static_simulation):
        result = static_simulation.stress_principal_elemental(components=3)
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("S3")
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (8,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises(self, static_simulation):
        result = static_simulation.stress_eqv_von_mises()
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("S_eqv")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (64,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_elemental(self, static_simulation):
        stress_vm = static_simulation.stress_eqv_von_mises_elemental()
        assert len(stress_vm._fc) == 1
        assert stress_vm._fc.get_time_scoping().ids == [1]
        field = stress_vm._fc[0]
        op = static_simulation._model.operator("S_eqv")
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (8,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_nodal(self, static_simulation):
        stress_vm = static_simulation.stress_eqv_von_mises_nodal()
        assert len(stress_vm._fc) == 1
        assert stress_vm._fc.get_time_scoping().ids == [1]
        field = stress_vm._fc[0]
        op = static_simulation._model.operator("S_eqv")
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (81,)
        assert np.allclose(field.data, field_ref.data)

    def test_reaction_force(self, static_simulation):
        reaction_force = static_simulation.reaction_force()
        assert len(reaction_force._fc) == 1
        assert reaction_force._fc.get_time_scoping().ids == [1]
        field = reaction_force._fc[0]
        op = static_simulation._model.operator("RF")
        field_ref = op.eval()[0]
        assert field.component_count == 3
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_2:
            assert field.data.shape == (21, 3)
        else:
            assert field.data.shape == (81, 3)
        assert np.allclose(field.data, field_ref.data)

    def test_elemental_volume(self, static_simulation):
        elemental_volume = static_simulation.elemental_volume()
        assert len(elemental_volume._fc) == 1
        assert elemental_volume._fc.get_time_scoping().ids == [1]
        field = elemental_volume._fc[0]
        op = static_simulation._model.operator("ENG_VOL")
        field_ref = op.eval()[0]
        # print(field_ref)
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_stiffness_matrix_energy(self, static_simulation):
        stiffness_matrix_energy = static_simulation.stiffness_matrix_energy()
        assert len(stiffness_matrix_energy._fc) == 1
        assert stiffness_matrix_energy._fc.get_time_scoping().ids == [1]
        field = stiffness_matrix_energy._fc[0]
        op = static_simulation._model.operator("ENG_SE")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_artificial_hourglass_energy(self, static_simulation):
        artificial_hourglass_energy = static_simulation.artificial_hourglass_energy()
        assert len(artificial_hourglass_energy._fc) == 1
        assert artificial_hourglass_energy._fc.get_time_scoping().ids == [1]
        field = artificial_hourglass_energy._fc[0]
        op = static_simulation._model.operator("ENG_AHO")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_thermal_dissipation_energy(self, static_simulation):
        thermal_dissipation_energy = static_simulation.thermal_dissipation_energy()
        assert len(thermal_dissipation_energy._fc) == 1
        assert thermal_dissipation_energy._fc.get_time_scoping().ids == [1]
        field = thermal_dissipation_energy._fc[0]
        op = static_simulation._model.operator("ENG_TH")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_kinetic_energy(self, static_simulation):
        kinetic_energy = static_simulation.kinetic_energy()
        assert len(kinetic_energy._fc) == 1
        assert kinetic_energy._fc.get_time_scoping().ids == [1]
        field = kinetic_energy._fc[0]
        op = static_simulation._model.operator("ENG_KE")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature(self, static_simulation):
        structural_temperature = static_simulation.structural_temperature()
        assert len(structural_temperature._fc) == 1
        assert structural_temperature._fc.get_time_scoping().ids == [1]
        field = structural_temperature._fc[0]
        op = static_simulation._model.operator("BFE")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (192,)
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature_nodal(self, static_simulation):
        structural_temperature_nodal = static_simulation.structural_temperature_nodal()
        assert len(structural_temperature_nodal._fc) == 1
        assert structural_temperature_nodal._fc.get_time_scoping().ids == [1]
        field = structural_temperature_nodal._fc[0]
        op = static_simulation._model.operator("BFE")
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (81,)
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature_elemental(self, static_simulation):
        structural_temperature_elemental = (
            static_simulation.structural_temperature_elemental()
        )
        assert len(structural_temperature_elemental._fc) == 1
        assert structural_temperature_elemental._fc.get_time_scoping().ids == [1]
        field = structural_temperature_elemental._fc[0]
        op = static_simulation._model.operator("BFE")
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_thermal_strain(self, allkindofcomplexity):
        static_simulation = post.StaticMechanicalSimulation(allkindofcomplexity)
        # thermal_strain
        result = static_simulation.thermal_strain(components=1)
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("ETHX")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (40016,)
        assert np.allclose(field.data, field_ref.data)
        # thermal_strain_eqv
        result = static_simulation.thermal_strain_eqv()
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("ETH_EQV")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (40016,)
        assert np.allclose(field.data, field_ref.data)
        # thermal_strain_principal
        result = static_simulation.thermal_strain_principal(components=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("ETH1")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (40016,)
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises(self, static_simulation):
        result = static_simulation.elastic_strain_eqv_von_mises(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=static_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = static_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        field_ref = equivalent_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_nodal(self, static_simulation):
        result = static_simulation.elastic_strain_eqv_von_mises_nodal(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=static_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = static_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = static_simulation._model.operator(name="to_nodal_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_elemental(self, static_simulation):
        result = static_simulation.elastic_strain_eqv_von_mises_elemental(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = static_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=static_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = static_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = static_simulation._model.operator(name="to_elemental_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_plot_stress_eqv_von_mises(self, static_simulation, tmp_path):
        result = static_simulation.stress_eqv_von_mises_nodal()
        result.plot()
        d = tmp_path / "stress_eqv"
        d.mkdir()
        result.plot(screenshot=d / "stress.png")
        os.path.exists(d / "stress.png")

    def test_external_layer(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.displacement(external_layer=True)
        assert len(result.index.mesh_index) == 81
        assert np.allclose(
            result.max(axis="node_ids").array,
            [2.76941713e-09, 2.76940199e-09, 4.10914311e-10],
        )
        result = static_simulation.displacement(set_ids=[1], external_layer=[1, 2, 3])
        assert len(result.index.mesh_index) == 44
        result = static_simulation.stress_principal_elemental(external_layer=[1, 2, 3])
        assert len(result.index.mesh_index) == 3
        result = static_simulation.elastic_strain_eqv_von_mises_elemental(
            external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        result = static_simulation.stress_principal_nodal(external_layer=[1, 2, 3])
        assert len(result.index.mesh_index) == 44
        result = static_simulation.elastic_strain_eqv_von_mises_nodal(
            external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 44

    def test_skin_layer(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.displacement(skin=True)
        assert len(result.index.mesh_index) == 74
        assert np.allclose(
            result.max(axis="node_ids").array,
            [2.76941713e-09, 2.76940199e-09, 4.10914311e-10],
        )

    def test_skin_layer2(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.displacement(set_ids=[1], skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 44

    def test_skin_layer3(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.elastic_strain_eqv_von_mises_elemental(
            skin=[1, 2, 3]
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 14
        else:
            assert len(result.index.mesh_index) == 18

    def test_skin_layer4(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.stress_principal_nodal(skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 44

    def test_skin_layer5(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.elastic_strain_eqv_von_mises_nodal(skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 44

    def test_skin_layer6(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.stress_principal_elemental(skin=[1, 2, 3])
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 14
        else:
            assert len(result.index.mesh_index) == 18


# List of element configurations for each simulation type
element_configurations = {
    "static_simulation": {
        1: [1],
        2: [1, 2],
        3: [1, 2, 3],
        4: [1, 2, 3, 4],
        5: [1, 2, 3, 4, 5],
        6: [1, 2, 3, 4, 5, 6, 7, 8],
    },
    "transient_simulation": {
        1: [1],
        2: [1, 2],
        3: [1, 2, 3],
        4: [1, 2, 3, 4],
        5: [1, 2, 3, 4, 5],
        6: [1, 2, 3, 4, 5, 6, 7, 8],
    },
    "modal_simulation": {1: [1], 2: [1, 2], 3: [1, 2, 3], 4: [1, 2, 3, 4, 5, 6, 7, 8]},
    "harmonic_simulation": {1: [1], 2: [1, 2], 3: [1, 2, 3], 4: list(range(1, 100))},
    "cyclic_static_simulation": {
        # Empty dict because element selection is
        # not supported for cyclic simulations
    },
}

# Get a set of all element configurations defined in the dictionary above
all_configuration_ids = [True] + list(
    set().union(
        *[
            element_configurations.keys()
            for element_configurations in element_configurations.values()
        ]
    )
)


def compute_number_of_expected_nodes(on_skin: bool, average_per_body: bool):
    n_nodes_per_side = 4
    if on_skin:
        # Take all the surfaces and remove nodes at the edges
        # and corners that are counted 2 or 3 times.
        # Remove the edge that touches both bodies. It is counted
        # 3 times
        nodes_all_surfaces = n_nodes_per_side**2 * 7
        duplicate_nodes_on_edges = 11 * (n_nodes_per_side - 2)
        triplicated_nodes_at_corners = 7
        expected_number_of_nodes = (
            nodes_all_surfaces
            - duplicate_nodes_on_edges
            - 2 * triplicated_nodes_at_corners
            - 2 * n_nodes_per_side
        )
    else:
        n_solid_nodes = n_nodes_per_side**3
        n_shell_nodes_without_touching = (n_nodes_per_side - 1) * n_nodes_per_side
        expected_number_of_nodes = n_solid_nodes + n_shell_nodes_without_touching

    if average_per_body:
        # Add boundary nodes again (duplicate nodes at the boundary)
        expected_number_of_nodes += n_nodes_per_side

    return expected_number_of_nodes


def get_shell_scoping(solid_mesh: MeshedRegion):
    split_scoping = operators.scoping.split_on_property_type()
    split_scoping.inputs.mesh(solid_mesh)
    split_scoping.inputs.label1("mat")
    split_scoping.inputs.requested_location(locations.elemental)

    splitted_scoping = split_scoping.eval()

    return splitted_scoping.get_scoping({"mat": 2})


def _check_nodal_across_body_results(
    fields_container: FieldsContainer,
    expected_results: dict[int, dict[str, float]],
    on_skin: bool,
):
    number_of_nodes_checked = 0
    assert len(fields_container) == 1
    field = fields_container[0]
    for node_id, expected_result_per_node in expected_results.items():
        if node_id in field.scoping.ids:
            number_of_nodes_checked += 1
            actual_result = field.get_entity_data_by_id(node_id)

            values_for_node = np.array(list(expected_result_per_node.values()))
            assert values_for_node.size > 0
            assert values_for_node.size < 3
            avg_expected_result = np.mean(values_for_node)

            if on_skin and len(values_for_node) > 1:
                # Skip elements at the edge that connects the body
                # because the averaging on the skin is different. For instance
                # 3 skin elements are involved the averaging of the inner elements
                continue
            assert np.isclose(
                actual_result, avg_expected_result, rtol=1e-3
            ), f"{values_for_node}, {node_id}"
    return number_of_nodes_checked


def _check_nodal_average_per_body_results(
    fields_container: FieldsContainer,
    expected_results: dict[int, dict[str, float]],
):
    number_of_nodes_checked = 0
    for node_id, expected_result_per_node in expected_results.items():
        for material in [1, 2]:
            field = fields_container.get_field({"mat": material})
            if node_id in field.scoping.ids:
                number_of_nodes_checked += 1
                actual_result = field.get_entity_data_by_id(node_id)
                expected_result = expected_result_per_node[str(material)]
                assert np.isclose(actual_result, expected_result, rtol=1e-3)
    return number_of_nodes_checked


def _check_elemental_per_body_results(
    fields_container: FieldsContainer,
    expected_results: dict[int, float],
    shell_elements_scoping: Scoping,
    element_id_to_skin_ids: dict[int, list[int]],
):
    checked_elements = 0

    for element_id, expected_value in expected_results.items():
        if element_id not in shell_elements_scoping.ids:
            continue
        skin_ids = element_id_to_skin_ids[element_id]
        for skin_id in skin_ids:
            for material in [1, 2]:
                field = fields_container.get_field({"mat": material})
                if skin_id in field.scoping.ids:
                    assert np.isclose(
                        field.get_entity_data_by_id(skin_id),
                        expected_value,
                        rtol=1e-3,
                    )
                    checked_elements += 1
    return checked_elements


def _check_elemental_across_body_results(
    fields_container: FieldsContainer,
    expected_results: dict[int, float],
    shell_elements_scoping: Scoping,
    element_id_to_skin_ids: dict[int, list[int]],
):
    checked_elements = 0

    for element_id, expected_value in expected_results.items():
        if element_id not in shell_elements_scoping.ids:
            continue
        skin_ids = element_id_to_skin_ids[element_id]
        for skin_id in skin_ids:
            assert len(fields_container) == 1
            field = fields_container[0]
            if skin_id in field.scoping.ids:
                assert np.isclose(
                    field.get_entity_data_by_id(skin_id),
                    expected_value,
                    rtol=1e-3,
                )
                checked_elements += 1
    return checked_elements


def _get_element_id_to_skin_id_map(skin_mesh: MeshedRegion, solid_mesh: MeshedRegion):
    skin_to_element_indices = skin_mesh.property_field("facets_to_ele")

    element_id_to_skin_ids = {}
    for skin_id in skin_mesh.elements.scoping.ids:
        element_idx = skin_to_element_indices.get_entity_data_by_id(skin_id)[0]
        solid_element_id = solid_mesh.elements.scoping.ids[element_idx]
        if solid_element_id not in element_id_to_skin_ids:
            element_id_to_skin_ids[solid_element_id] = []
        element_id_to_skin_ids[solid_element_id].append(skin_id)
    return element_id_to_skin_ids


# @pytest.mark.parametrize("average_per_body", [False, True])
# @pytest.mark.parametrize("on_skin", [True, False])
# # Note: shell_layer selection with multiple layers (e.g top/bottom)
# currently not working correctly
# # for mixed models.
# @pytest.mark.parametrize("shell_layer", [shell_layers.top, shell_layers.bottom])
# @pytest.mark.parametrize("location", [locations.elemental, locations.nodal])
# def test_shell_layer_extraction(
#     mixed_shell_solid_simulation,
#     shell_layer_multi_body_ref,
#     average_per_body,
#     on_skin,
#     shell_layer,
#     location,
# ):
#     if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1:
#         return
#
#     shell_layer_names = {shell_layers.top: "top", shell_layers.bottom: "bot"}
#
#     if average_per_body:
#         averaging_config = AveragingConfig(
#             body_defining_properties=["mat"], average_per_body=True
#         )
#     else:
#         averaging_config = AveragingConfig(
#             body_defining_properties=None, average_per_body=False
#         )
#
#     res = mixed_shell_solid_simulation._get_result(
#         base_name="S",
#         skin=on_skin,
#         components=["X"],
#         location=location,
#         category=ResultCategory.matrix,
#         shell_layer=shell_layer,
#         averaging_config=averaging_config,
#     )
#
#     if location == locations.nodal:
#         expected_results = get_ref_per_body_results_mechanical(
#             shell_layer_multi_body_ref[
#                 f"stress_{shell_layer_names[shell_layer]}_nodal"
#             ],
#             mixed_shell_solid_simulation.mesh._meshed_region,
#         )
#
#         expected_number_of_nodes = compute_number_of_expected_nodes(
#             on_skin, average_per_body
#         )
#
#         if average_per_body:
#             number_of_nodes_checked = _check_nodal_average_per_body_results(
#                 fields_container=res._fc,
#                 expected_results=expected_results,
#             )
#         else:
#             number_of_nodes_checked = _check_nodal_across_body_results(
#                 fields_container=res._fc,
#                 expected_results=expected_results,
#                 on_skin=on_skin,
#             )
#
#         assert number_of_nodes_checked == expected_number_of_nodes
#
#     else:
#         ref_result = get_ref_result_per_element(
#             shell_layer_multi_body_ref[
#                 f"stress_{shell_layer_names[shell_layer]}_elemental"
#             ].combined
#         )
#         checked_elements = 0
#
#         if on_skin:
#             skin_mesh = res._fc[0].meshed_region
#             solid_mesh = mixed_shell_solid_simulation.mesh._meshed_region
#
#             shell_elements_scoping = get_shell_scoping(solid_mesh)
#             element_id_to_skin_ids = _get_element_id_to_skin_id_map(
#                 skin_mesh, solid_mesh
#             )
#
#             # Note: In this branch only shell elements are checked,
#             # since only the shell elements are
#             # affected by the shell layer extraction.
#             # The skin of the solid elements is cumbersome to
#             # extract and check and is skipped here.
#             if average_per_body:
#                 checked_elements = _check_elemental_per_body_results(
#                     fields_container=res._fc,
#                     expected_results=ref_result,
#                     shell_elements_scoping=shell_elements_scoping,
#                     element_id_to_skin_ids=element_id_to_skin_ids,
#                 )
#             else:
#                 checked_elements = _check_elemental_across_body_results(
#                     fields_container=res._fc,
#                     expected_results=ref_result,
#                     shell_elements_scoping=shell_elements_scoping,
#                     element_id_to_skin_ids=element_id_to_skin_ids,
#                 )
#
#             assert checked_elements == 9
#         else:
#             for element_id, expected_value in ref_result.items():
#                 if average_per_body:
#                     for material in [1, 2]:
#                         field = res._fc.get_field({"mat": material})
#                         if element_id in field.scoping.ids:
#                             assert np.isclose(
#                                 field.get_entity_data_by_id(element_id),
#                                 expected_value,
#                                 rtol=1e-3,
#                             ), expected_value
#                             checked_elements += 1
#                 else:
#                     assert np.isclose(
#                         res._fc[0].get_entity_data_by_id(element_id),
#                         expected_value,
#                         rtol=1e-3,
#                     ), expected_value
#                     checked_elements += 1
#             assert checked_elements == 36


@pytest.mark.parametrize("average_per_body", [False, True])
@pytest.mark.parametrize("on_skin", [True, False])
# Note: shell_layer selection with multiple layers (e.g top/bottom) currently not working correctly
# for mixed models.
@pytest.mark.parametrize("shell_layer", [shell_layers.top, shell_layers.bottom])
@pytest.mark.parametrize("location", [locations.elemental, locations.nodal])
@pytest.mark.parametrize(
    "simulation_str",
    ["modal_simulation", "transient_simulation", "harmonic_simulation"],
)
def test_shell_layer_extraction_all_types_except_static(
    on_skin, location, shell_layer, average_per_body, simulation_str, request
):
    # This test just verifies that the shape of the results is correct for all simulation types,
    # to ensure the arguments are passed correctly.
    # The correctness of the results is not checked here, but for StaticMechanicalSimulation, since
    # the workflow is the same for all simulation types.
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1:
        return

    if average_per_body:
        averaging_config = AveragingConfig(
            body_defining_properties=["mat"], average_per_body=True
        )
    else:
        averaging_config = AveragingConfig(
            body_defining_properties=None, average_per_body=False
        )

    simulation = request.getfixturevalue(simulation_str)

    res = simulation._get_result(
        base_name="S",
        skin=on_skin,
        components=["X"],
        location=location,
        category=ResultCategory.matrix,
        shell_layer=shell_layer,
        averaging_config=averaging_config,
    )

    for field in res._fc:
        for entity_id in field.scoping.ids:
            entity_data = field.get_entity_data_by_id(entity_id)
            assert entity_data.shape == (1,)


@pytest.mark.parametrize(
    "average_per_body",
    [
        False,
        pytest.param(
            True,
            marks=pytest.mark.xfail(
                reason="Failing because scopings without results"
                " are not handled correctly in the current implementation."
            ),
        ),
    ],
)
@pytest.mark.parametrize("on_skin", [True, False])
# Note: shell_layer selection with multiple layers (e.g top/bottom) currently not working correctly
# for mixed models.
@pytest.mark.parametrize("shell_layer", [shell_layers.top, shell_layers.bottom])
@pytest.mark.parametrize("location", [locations.elemental, locations.nodal])
@pytest.mark.parametrize(
    "simulation_str",
    [
        "two_cubes_contact_simulation",
        pytest.param(
            "mixed_shell_solid_with_contact_simulation",
            marks=pytest.mark.xfail(
                reason="Failing because scopings without results"
                " are not handled correctly in the current implementation."
            ),
        ),
    ],
)
def test_shell_layer_extraction_contacts(
    simulation_str, average_per_body, on_skin, shell_layer, location, request
):
    # Test some models with contacts, because models with contacts
    # result in fields without results which can cause problems in conjunction
    # with shell layer extraction.
    simulation = request.getfixturevalue(simulation_str)

    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1:
        return

    if average_per_body:
        averaging_config = AveragingConfig(
            body_defining_properties=["mat"], average_per_body=True
        )
    else:
        averaging_config = AveragingConfig(
            body_defining_properties=None, average_per_body=False
        )

    res = simulation._get_result(
        base_name="S",
        skin=on_skin,
        components=["X"],
        location=location,
        category=ResultCategory.equivalent,
        shell_layer=shell_layer,
        averaging_config=averaging_config,
    )

    # Just do a rough comparison.
    # This test is mainly to check if the
    # workflow runs without errors because of
    # empty fields for some materials
    max_val = res.max().array[0]
    if simulation_str == "two_cubes_contact_simulation":
        assert max_val > 1 and max_val < 1.1
    else:
        assert max_val > 7.7 and max_val < 7.8


@pytest.mark.parametrize("skin", all_configuration_ids)
@pytest.mark.parametrize("result_name", ["stress", "elastic_strain", "displacement"])
@pytest.mark.parametrize("mode", [None, "principal", "equivalent"])
@pytest.mark.parametrize(
    "simulation_str",
    [
        "static_simulation",
        "transient_simulation",
        "modal_simulation",
        "harmonic_simulation",
        # Just some very basic tests for the cyclic simulation
        "cyclic_static_simulation",
    ],
)
def test_skin_extraction(skin, result_name, mode, simulation_str, request):
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1:
        return

    time_id = 1

    simulation = request.getfixturevalue(simulation_str)

    supports_elemental = True
    is_cyclic_simulation = simulation_str == "cyclic_static_simulation"

    if is_cyclic_simulation:
        if result_name == "elastic_strain":
            # cyclic simulation does not have elastic strain results
            return
        if is_equivalent(mode) or is_principal(mode):
            # Test for equivalent and principal modes not implemented
            return

    if skin is not True:
        skin = element_configurations[simulation_str].get(skin)
        if skin is None:
            # Return if a element configuration does
            # not exist for a given simulation type
            return

    if result_name == "displacement":
        supports_elemental = False
        if is_principal(mode) or is_equivalent(mode):
            # Return for unsupported results
            return

    if isinstance(skin, list):
        element_ids = skin
    else:
        if isinstance(simulation, post.ModalMechanicalSimulation):
            # The modal result contains different element types. Here
            # we just extract the solid elements
            solid_elements_mesh = simulation.split_mesh_by_properties(
                {elemental_properties.element_type: element_types.Hex20.value}
            )
            if isinstance(solid_elements_mesh, Meshes):
                element_ids = solid_elements_mesh[0].element_ids
            else:
                element_ids = solid_elements_mesh.element_ids
            skin = element_ids
        else:
            element_ids = simulation.mesh.element_ids

    scoping = None
    if isinstance(skin, list):
        scoping = Scoping(ids=element_ids, location="elemental")

    fc_elemental_nodal = None
    if supports_elemental:
        fc_elemental_nodal = get_elemental_nodal_results(
            simulation=simulation,
            result_name=result_name,
            scoping=scoping,
            expand_cyclic=is_cyclic_simulation,
            mode=mode,
        )

        get_and_check_elemental_skin_results(
            static_simulation=simulation,
            fc_elemental_nodal=fc_elemental_nodal,
            result_name=result_name,
            mode=mode,
            element_ids=element_ids,
            skin=skin,
            expand_cyclic=is_cyclic_simulation,
        )

    # Not all the simulation types have the expand_cyclic argument
    kwargs = {}
    if is_cyclic_simulation:
        kwargs["expand_cyclic"] = True

    # For displacements the nodal result
    # is just called displacement without
    # the "nodal" suffix
    nodal_suffix = "_nodal"
    if result_name == "displacement":
        nodal_suffix = ""

    result_skin_scoped_nodal = getattr(
        simulation, f"{result_name}{mode_suffix(mode)}{nodal_suffix}"
    )(set_ids=[time_id], skin=skin, **kwargs)
    nodal_skin_field = result_skin_scoped_nodal._fc[0]

    expected_nodal_values_field = get_expected_nodal_skin_results(
        simulation=simulation,
        result_name=result_name,
        mode=mode,
        skin_mesh=nodal_skin_field.meshed_region,
        elemental_nodal_results=fc_elemental_nodal,
        expand_cyclic=is_cyclic_simulation,
    )

    for node_id in expected_nodal_values_field.scoping.ids:
        if result_name == "displacement":
            if node_id not in nodal_skin_field.scoping.ids:
                # We get the displacement results also for internal
                # nodes. We skip these nodes here.
                continue
        assert np.allclose(
            expected_nodal_values_field.get_entity_data_by_id(node_id),
            nodal_skin_field.get_entity_data_by_id(node_id),
        ), str(node_id)

    # result_skin_scoped_elemental_nodal = getattr(
    #     static_simulation, f"{result_name}{mode_suffix(mode)}"
    # )(all_sets=True, skin=element_ids)

    # Todo: Elemental nodal does not work
    # Returns just the element nodal data of the solid
    # result_skin_scoped_elemental_nodal


class TestTransientMechanicalSimulation:
    def test_times_argument(self, transient_simulation, static_simulation):
        with pytest.raises(
            ValueError, match="Could not find time=0.0 in the simulation."
        ):
            _ = transient_simulation.displacement(times=0.0)

        # Get reference field at t=0.15s
        op = transient_simulation._model.operator("UX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            15, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        field_ref = op.eval()[0]
        # Test for times= exact float
        result = transient_simulation.displacement(components=["X"], times=0.15)
        field = result._fc[0]
        assert np.allclose(field.data, field_ref.data)
        # Test for times= near float
        result = transient_simulation.displacement(components=["X"], times=0.1496)
        field = result._fc[0]
        assert np.allclose(field.data, field_ref.data)
        # Test for times= just not near float
        with pytest.raises(
            ValueError, match="Could not find time=0.1495 in the simulation."
        ):
            _ = transient_simulation.displacement(components=["X"], times=0.1495)

    def test_displacement(self, transient_simulation):
        result = transient_simulation.displacement(
            components=["X"],
            node_ids=[2, 3, 4],
            all_sets=True,
        )
        assert len(result._fc) == 20
        assert len(result._fc.get_time_scoping().ids) == 20
        result = transient_simulation.displacement(components=["X"], node_ids=[2, 3, 4])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [20]
        field = result._fc[0]
        op = transient_simulation._model.operator("UX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            20, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [2, 3, 4], server=transient_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

        result = transient_simulation.displacement(
            components=1,
            named_selections=transient_simulation.named_selections[:2],
            set_ids=[2],
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        assert field.component_count == 1
        assert field.data.shape == (393,)

    def test_velocity(self, transient_simulation):
        result = transient_simulation.velocity(
            components=["X"], node_ids=[2, 3, 4], set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("VX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [2, 3, 4], server=transient_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

    def test_acceleration(self, transient_simulation):
        result = transient_simulation.acceleration(
            components=["X"], node_ids=[2, 3, 4], set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("AX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [2, 3, 4], server=transient_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress(self, transient_simulation):
        result = transient_simulation.stress(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_elemental(self, transient_simulation):
        result = transient_simulation.stress_elemental(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_nodal(self, transient_simulation):
        result = transient_simulation.stress_nodal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal(self, transient_simulation):
        result = transient_simulation.stress_principal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S1")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_nodal(self, transient_simulation):
        result = transient_simulation.stress_principal_nodal(components=2, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S2")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_elemental(self, transient_simulation):
        result = transient_simulation.stress_principal_elemental(
            components=3, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S3")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises(self, transient_simulation):
        result = transient_simulation.stress_eqv_von_mises(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_elemental(self, transient_simulation):
        result = transient_simulation.stress_eqv_von_mises_elemental(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_nodal(self, transient_simulation):
        result = transient_simulation.stress_eqv_von_mises_nodal(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain(self, transient_simulation):
        result = transient_simulation.elastic_strain(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_elemental(self, transient_simulation):
        result = transient_simulation.elastic_strain_elemental(
            components=1, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_nodal(self, transient_simulation):
        result = transient_simulation.elastic_strain_nodal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal(self, transient_simulation):
        result = transient_simulation.elastic_strain_principal(
            components=1, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        principal_op = transient_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_1()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_nodal(self, transient_simulation):
        result = transient_simulation.elastic_strain_principal_nodal(
            components=2, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        principal_op = transient_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_2()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_elemental(self, transient_simulation):
        result = transient_simulation.elastic_strain_principal_elemental(
            components=3, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        principal_op = transient_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_3()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_reaction_force(self, allkindofcomplexity):
        transient_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.reaction_force(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("RF")
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    def test_elemental_volume(self, transient_simulation):
        result = transient_simulation.elemental_volume(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENG_VOL")
        field_ref = op.eval()[0]
        # print(field_ref)
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_artificial_hourglass_energy(self, transient_simulation):
        result = transient_simulation.artificial_hourglass_energy(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENG_AHO")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_thermal_dissipation_energy(self, transient_simulation):
        result = transient_simulation.thermal_dissipation_energy(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENG_TH")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_thermal_strain(self, allkindofcomplexity):
        simulation = post.TransientMechanicalSimulation(allkindofcomplexity)
        print(simulation)
        # thermal_strain
        result = simulation.thermal_strain(components=1)
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = simulation._model.operator("ETHX")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (40016,)
        assert np.allclose(field.data, field_ref.data)
        # thermal_strain_eqv
        result = simulation.thermal_strain_eqv()
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = simulation._model.operator("ETH_EQV")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (40016,)
        assert np.allclose(field.data, field_ref.data)
        # thermal_strain_principal
        result = simulation.thermal_strain_principal(components=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = simulation._model.operator("ETH1")
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (40016,)
        assert np.allclose(field.data, field_ref.data)

    def test_kinetic_energy(self, transient_simulation):
        result = transient_simulation.kinetic_energy(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENG_KE")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature(self, transient_simulation):
        result = transient_simulation.structural_temperature(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("BFE")
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature_nodal(self, transient_simulation):
        result = transient_simulation.structural_temperature_nodal(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("BFE")
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature_elemental(self, transient_simulation):
        result = transient_simulation.structural_temperature_elemental(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("BFE")
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces(self, allkindofcomplexity):
    #     transient_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.transient_mechanical,
    #     )
    #     result = transient_simulation.element_nodal_forces(set_ids=[1])
    #     assert len(result._fc) == 1
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = transient_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)
    #
    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces_nodal(self, allkindofcomplexity):
    #     transient_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.transient_mechanical,
    #     )
    #     result = transient_simulation.element_nodal_forces_nodal(set_ids=[1])
    #     assert len(result._fc) == 3
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = transient_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     op.connect(9, post.locations.nodal)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)
    #
    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces_elemental(self, allkindofcomplexity):
    #     transient_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.transient_mechanical,
    #     )
    #     result = transient_simulation.element_nodal_forces_elemental(set_ids=[1])
    #     assert len(result._fc) == 3
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = transient_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     op.connect(9, post.locations.elemental)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises(self, transient_simulation):
        result = transient_simulation.elastic_strain_eqv_von_mises(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = transient_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        field_ref = equivalent_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_nodal(self, transient_simulation):
        result = transient_simulation.elastic_strain_eqv_von_mises_nodal(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = transient_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = transient_simulation._model.operator(name="to_nodal_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_elemental(self, transient_simulation):
        result = transient_simulation.elastic_strain_eqv_von_mises_elemental(
            set_ids=[1]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = transient_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = transient_simulation._model.operator(name="to_elemental_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_external_layer(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.displacement(all_sets=True, external_layer=True)
        assert len(result.columns.set_ids) == 20
        assert len(result.index.mesh_index) == 393
        assert np.allclose(
            result.select(set_ids=[2]).max(axis="node_ids").array,
            [5.14806800e-07, 1.63151192e-03, 9.78100326e-06],
        )
        result = transient_simulation.displacement(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 44
        result = transient_simulation.stress_principal_elemental(
            external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        result = transient_simulation.elastic_strain_eqv_von_mises_elemental(
            external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        result = transient_simulation.stress_principal_nodal(external_layer=[1, 2, 3])
        assert len(result.index.mesh_index) == 44
        result = transient_simulation.elastic_strain_eqv_von_mises_nodal(
            external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 44

    def test_skin_layer(self, transient_simulation: post.TransientMechanicalSimulation):
        result = transient_simulation.displacement(all_sets=True, skin=True)
        assert len(result.columns.set_ids) == 20
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 374
        else:
            assert len(result.index.mesh_index) == 393
        assert np.allclose(
            result.select(set_ids=[2]).max(axis="node_ids").array,
            [5.14806800e-07, 1.63151192e-03, 9.78100326e-06],
        )

    def test_skin_layer2(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.displacement(set_ids=[1], skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 44

    def test_skin_layer3(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.stress_principal_elemental(
            skin=list(range(1, 100))
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 124
        else:
            assert len(result.index.mesh_index) == 240

    def test_skin_layer4(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.elastic_strain_eqv_von_mises_elemental(
            skin=list(range(1, 100))
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 124
        else:
            assert len(result.index.mesh_index) == 240

    def test_skin_layer5(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.stress_principal_nodal(skin=list(range(1, 100)))
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 374
        else:
            assert len(result.index.mesh_index) == 393

    def test_skin_layer6(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.elastic_strain_eqv_von_mises_nodal(
            skin=list(range(1, 100))
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 374
        else:
            assert len(result.index.mesh_index) == 393


class TestModalMechanicalSimulation:
    @fixture
    def frame_modal_simulation(self, modalframe):
        return post.load_simulation(
            data_sources=modalframe,
            simulation_type=AvailableSimulationTypes.modal_mechanical,
        )

    def test_cyclic(self, simple_cyclic):
        simulation = post.ModalMechanicalSimulation(simple_cyclic)
        displacement = simulation.displacement(expand_cyclic=False)
        assert "base_sector" in displacement.columns.names
        assert len(displacement.mesh_index) == 51

        displacement = simulation.displacement(expand_cyclic=True)
        assert "base_sector" not in displacement.columns.names
        assert len(displacement.mesh_index) == 408

        with pytest.raises(
            ValueError,
            match="'phase_angle_cyclic' argument only accepts a single float value.",
        ):
            _ = simulation.displacement(phase_angle_cyclic=[0.1])

        with pytest.raises(
            ValueError,
            match="Sector selection with 'expand_cyclic' starts at 1.",
        ):
            _ = simulation.displacement(expand_cyclic=[0])

        displacement = simulation.displacement(phase_angle_cyclic=90)
        assert displacement
        displacement = simulation.displacement(phase_angle_cyclic=90.0)
        assert displacement

    def test_multi_stage(self, multi_stage_cyclic):
        simulation = post.ModalMechanicalSimulation(multi_stage_cyclic)

        displacement = simulation.displacement(expand_cyclic=False)
        assert "base_sector" in displacement.columns.names
        assert "stage" in displacement.columns.names
        assert len(displacement.mesh_index) == 3595

        displacement = simulation.displacement(expand_cyclic=True)
        assert "base_sector" not in displacement.columns.names
        assert "stage" not in displacement.columns.names
        assert len(displacement.mesh_index) == 26742

        with pytest.raises(
            ValueError,
            match="Sector selection with 'expand_cyclic' starts at 1.",
        ):
            _ = simulation.displacement(expand_cyclic=[[0, 1], 1])

        displacement = simulation.displacement(expand_cyclic=[1, 2])
        assert "base_sector" not in displacement.columns.names
        assert "stage" not in displacement.columns.names
        assert len(displacement.mesh_index) == 18717

        displacement = simulation.displacement(expand_cyclic=[[1, 2], 1])
        assert "base_sector" not in displacement.columns.names
        assert "stage" not in displacement.columns.names
        assert len(displacement.mesh_index) == 5644

        displacement = simulation.displacement(expand_cyclic=[[1, 2], [1, 2]])
        # print(displacement)
        assert "base_sector" not in displacement.columns.names
        assert "stage" not in displacement.columns.names
        assert len(displacement.mesh_index) == 6848

        with pytest.raises(
            ValueError, match="'expand_cyclic' only accepts lists of int values >= 1."
        ):
            _ = simulation.displacement(expand_cyclic=[[1, 2], [0.2, 2]])

        with pytest.raises(
            ValueError,
            match="'expand_cyclic' argument can only be a boolean or a list.",
        ):
            _ = simulation.displacement(expand_cyclic=1)

    def test_displacement(self, modal_simulation):
        # print(modal_simulation)
        result = modal_simulation.displacement(
            components=["X"],
            node_ids=[2, 3, 4],
            all_sets=True,
        )
        assert len(result._fc) == 45
        assert len(result._fc.get_time_scoping().ids) == 45

        result = modal_simulation.displacement(components=["X"], node_ids=[2, 3, 4])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("UX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [2, 3, 4], server=modal_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

    def test_reaction_force(self, allkindofcomplexity):
        modal_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.modal_mechanical,
        )
        result = modal_simulation.reaction_force(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("RF")
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces(self, allkindofcomplexity):
    #     modal_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.modal_mechanical,
    #     )
    #     result = modal_simulation.element_nodal_forces(set_ids=[1])
    #     assert len(result._fc) == 1
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = modal_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)
    #
    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces_nodal(self, allkindofcomplexity):
    #     modal_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.modal_mechanical,
    #     )
    #     result = modal_simulation.element_nodal_forces_nodal(set_ids=[1])
    #     assert len(result._fc) == 3
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = modal_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     op.connect(9, post.locations.nodal)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)
    #
    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces_elemental(self, allkindofcomplexity):
    #     modal_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.modal_mechanical,
    #     )
    #     result = modal_simulation.element_nodal_forces_elemental(set_ids=[1])
    #     assert len(result._fc) == 3
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = modal_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     op.connect(9, post.locations.elemental)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)

    def test_stress(self, modal_simulation):
        result = modal_simulation.stress(components=1, modes=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_elemental(self, modal_simulation):
        result = modal_simulation.stress_elemental(components=1, set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_nodal(self, modal_simulation):
        result = modal_simulation.stress_nodal(components=1, set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal(self, modal_simulation):
        result = modal_simulation.stress_principal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("S1")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_nodal(self, modal_simulation):
        result = modal_simulation.stress_principal_nodal(components=2, set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("S2")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_elemental(self, modal_simulation):
        result = modal_simulation.stress_principal_elemental(components=3, set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("S3")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises(self, modal_simulation):
        result = modal_simulation.stress_eqv_von_mises(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_elemental(self, modal_simulation):
        result = modal_simulation.stress_eqv_von_mises_elemental(set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_nodal(self, modal_simulation):
        result = modal_simulation.stress_eqv_von_mises_nodal(set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elemental_volume(self, modal_simulation):
        result = modal_simulation.elemental_volume(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("ENG_VOL")
        field_ref = op.eval()[0]
        # print(field_ref)
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain(self, modal_simulation):
        result = modal_simulation.elastic_strain(components=1, modes=2)
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_elemental(self, modal_simulation):
        result = modal_simulation.elastic_strain_elemental(components=1, set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_nodal(self, modal_simulation):
        result = modal_simulation.elastic_strain_nodal(components=1, set_ids=[2])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal(self, modal_simulation):
        result = modal_simulation.elastic_strain_principal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        principal_op = modal_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_1()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_nodal(self, modal_simulation):
        result = modal_simulation.elastic_strain_principal_nodal(
            components=2, set_ids=[2]
        )
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        principal_op = modal_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_2()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_elemental(self, modal_simulation):
        result = modal_simulation.elastic_strain_principal_elemental(
            components=3, set_ids=[2]
        )
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            2, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        principal_op = modal_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_3()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises(self, modal_simulation):
        result = modal_simulation.elastic_strain_eqv_von_mises(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = modal_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        field_ref = equivalent_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_nodal(self, modal_simulation):
        result = modal_simulation.elastic_strain_eqv_von_mises_nodal(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = modal_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = modal_simulation._model.operator(name="to_nodal_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_elemental(self, modal_simulation):
        result = modal_simulation.elastic_strain_eqv_von_mises_elemental(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=modal_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = modal_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = modal_simulation._model.operator(name="to_elemental_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_disp_external_layer(
        self, frame_modal_simulation: post.ModalMechanicalSimulation
    ):
        result = frame_modal_simulation.displacement(set_ids=[1], external_layer=True)
        result_all = frame_modal_simulation.displacement(
            set_ids=[1], external_layer=False
        )
        assert len(result.index.mesh_index) == 5886
        assert len(result.index.mesh_index) == len(result_all.index.mesh_index)
        assert np.allclose(
            result.max(axis="node_ids").array, [0.05656421, 9.59989137, 1.08656671]
        )
        result = frame_modal_simulation.displacement(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 21
        assert np.allclose(
            result.max(axis="node_ids").array, [-0.77876072, 7.08211902, 0.05292333]
        )

    def test_stress_external_layer(
        self, frame_modal_simulation: post.ModalMechanicalSimulation
    ):
        result = frame_modal_simulation.stress_elemental(
            all_sets=True, external_layer=True
        )
        assert len(result.index.mesh_index) == 2842
        assert len(result.columns.set_ids) == 6
        assert np.allclose(
            result.select(set_ids=[3]).max(axis="element_ids").array,
            [
                [
                    464.27737236,
                    627.19576979,
                    1661.52572632,
                    285.47153473,
                    682.4336586,
                    501.27880096,
                ]
            ],
        )
        result = frame_modal_simulation.stress_elemental(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        assert len(result.columns.set_ids) == 1
        assert np.allclose(
            result.max(axis="element_ids").array,
            [
                [
                    -6.46462837,
                    1.80925381,
                    107.2106514,
                    2.10453892,
                    9.44412744,
                    -2.84251213,
                ]
            ],
        )
        result = frame_modal_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 21
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [280.45842315])

    def test_stress_external_layer2(
        self, frame_modal_simulation: post.ModalMechanicalSimulation
    ):
        result = frame_modal_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], external_layer=True
        )
        assert len(result.index.mesh_index) == 5886
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [1285.17926915])

    def test_strain_external_layer(
        self, frame_modal_simulation: post.ModalMechanicalSimulation
    ):
        result = frame_modal_simulation.stress_principal_elemental(
            all_sets=True, external_layer=True
        )
        assert len(result.index.mesh_index) == 2842
        assert len(result.columns.set_ids) == 6
        assert np.allclose(
            result.select(set_ids=[1]).max(axis="element_ids").array, [1282.65478454]
        )
        result = frame_modal_simulation.stress_principal_elemental(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="element_ids").array, [123.69229739])
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 21
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [0.00149252])

    def test_strain_external_layer2(
        self, frame_modal_simulation: post.ModalMechanicalSimulation
    ):
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], external_layer=True
        )
        assert len(result.index.mesh_index) == 5886
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [0.00684776])
        result = frame_modal_simulation.elastic_strain_principal_nodal(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 21
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [0.00131161])

    def test_strain_external_layer3(
        self, frame_modal_simulation: post.ModalMechanicalSimulation
    ):
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_elemental(
            set_ids=[1], external_layer=True
        )
        assert len(result.index.mesh_index) == 2842
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="element_ids").array, [0.00620539])

    def test_disp_skin(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        result = frame_modal_simulation.displacement(set_ids=[1], skin=True)
        result_all = frame_modal_simulation.displacement(set_ids=[1], skin=False)
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 4068
        else:
            assert len(result.index.mesh_index) == 5828
        assert np.allclose(
            result.max(axis="node_ids").array, [0.05656421, 9.59989137, 1.08656671]
        )
        result = frame_modal_simulation.displacement(set_ids=[1], skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 21
        assert np.allclose(
            result.max(axis="node_ids").array, [-0.77876072, 7.08211902, 0.05292333]
        )

    def test_stress_skin(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        if frame_modal_simulation._model._server.meet_version("7.1"):
            result = frame_modal_simulation.stress_elemental(all_sets=True, skin=True)
            assert len(result.index.mesh_index) == 2048
            assert len(result.columns.set_ids) == 6
        elif frame_modal_simulation._model._server.meet_version("6.2"):
            result = frame_modal_simulation.stress_elemental(all_sets=True, skin=True)
            assert len(result.index.mesh_index) == 11146
            assert len(result.columns.set_ids) == 6
        result = frame_modal_simulation.stress_elemental(
            set_ids=[1], skin=list(range(1, 100))
        )
        assert len(result.columns.set_ids) == 1
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
            assert len(result.index.mesh_index) == 132
            assert np.allclose(
                result.max(axis="element_ids").array,
                [
                    [
                        88.09000492095947,
                        426.211181640625,
                        747.8219401041666,
                        30.50066868464152,
                        412.8089192708333,
                        109.25983428955078,
                    ]
                ],
            )
        elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 36
            assert np.allclose(
                result.max(axis="element_ids").array,
                [
                    [
                        36.52192259,
                        58.73246002,
                        371.72294617,
                        12.80614456,
                        134.60557556,
                        38.0447108,
                    ]
                ],
            )
        else:
            assert len(result.index.mesh_index) == 110
            assert np.allclose(
                result.max(axis="element_ids").array,
                [
                    [
                        36.52192259,
                        58.73246002,
                        371.72294617,
                        25.97949378,
                        139.83338165,
                        69.25232569,
                    ]
                ],
            )

    def test_stress_skin2(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        result = frame_modal_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], skin=frame_modal_simulation.mesh.element_ids
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 4068
            assert np.allclose(result.max(axis="node_ids").array, [1295.83764693])
        else:
            assert len(result.index.mesh_index) == 5828
            assert np.allclose(result.max(axis="node_ids").array, [1285.17926915])
        assert len(result.columns.set_ids) == 1
        result = frame_modal_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], skin=True
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 4068
        else:
            assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1

    def test_strain_skin(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        if frame_modal_simulation._model._server.meet_version("7.1"):
            result = frame_modal_simulation.stress_principal_elemental(
                all_sets=True, skin=True
            )
            assert len(result.index.mesh_index) == 2048
            assert len(result.columns.set_ids) == 6
            assert np.allclose(
                result.select(set_ids=[1]).max(axis="element_ids").array,
                [1339.75343629],
            )
        elif frame_modal_simulation._model._server.meet_version("6.2"):
            result = frame_modal_simulation.stress_principal_elemental(
                all_sets=True, skin=True
            )
            assert len(result.index.mesh_index) == 11146
            assert len(result.columns.set_ids) == 6
            assert np.allclose(
                result.select(set_ids=[1]).max(axis="element_ids").array,
                [1602.16293782],
            )
        result = frame_modal_simulation.stress_principal_elemental(
            set_ids=[1], skin=list(range(1, 100))
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
            assert len(result.index.mesh_index) == 132
        elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 36
        else:
            assert len(result.index.mesh_index) == 110
        assert len(result.columns.set_ids) == 1

    def test_strain_skin2(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], skin=frame_modal_simulation.mesh.element_ids
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 4068
            assert np.allclose(result.max(axis="node_ids").array, [0.00695066])
        else:
            assert len(result.index.mesh_index) == 5828
            assert np.allclose(result.max(axis="node_ids").array, [0.00684776])
        assert len(result.columns.set_ids) == 1

    def test_strain_skin3(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], skin=True
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 4068
        else:
            assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1
        result = frame_modal_simulation.elastic_strain_principal_nodal(
            set_ids=[1], skin=frame_modal_simulation.mesh.element_ids
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 4068
        else:
            assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_elemental(
            set_ids=[1], skin=True
        )
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
            assert len(result.index.mesh_index) == 2048
        else:
            assert len(result.index.mesh_index) == 11146
        assert len(result.columns.set_ids) == 1


class TestHarmonicMechanicalSimulation:
    def test_cyclic(self, simple_cyclic):
        simulation = post.HarmonicMechanicalSimulation(simple_cyclic)
        result = simulation.displacement(expand_cyclic=False)
        # print(result)
        assert "base_sector" in result.columns.names
        result = simulation.displacement(expand_cyclic=True)
        # print(result)
        assert "base_sector" not in result.columns.names

    def test_multi_stage(self, multi_stage_cyclic):
        simulation = post.HarmonicMechanicalSimulation(multi_stage_cyclic)
        result = simulation.displacement(expand_cyclic=False)
        # print(result)
        assert "base_sector" in result.columns.names
        assert "stage" in result.columns.names
        result = simulation.displacement(expand_cyclic=True)
        # print(result)
        assert "base_sector" not in result.columns.names
        assert "stage" not in result.columns.names

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
        reason="Available starting DPF 4.0",
    )
    def test_with_grpc_server(self, complex_model, grpc_server):
        simulation = post.HarmonicMechanicalSimulation(
            complex_model, server=grpc_server
        )
        assert simulation._model._server != dpf.SERVER
        _ = simulation.displacement()
        _ = simulation.displacement(skin=True)

    def test_displacement(self, harmonic_simulation):
        # print(harmonic_simulation)

        result = harmonic_simulation.displacement(components=["X"], node_ids=[2, 3, 4])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("UX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
            [2, 3, 4], server=harmonic_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

    # def test_amplitude(self, harmonic_simulation):
    #     result = harmonic_simulation.displacement(
    #         components=["X"], node_ids=[2, 3, 4], amplitude=True
    #     )
    #     assert len(result._fc) == 1
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #
    #     op = harmonic_simulation._model.operator("UX")
    #     time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
    #         1, server=harmonic_simulation._model._server
    #     )
    #     op.connect(0, time_scoping)
    #     mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
    #         [2, 3, 4], server=harmonic_simulation._model._server
    #     )
    #     op.connect(1, mesh_scoping)
    #     amplitude_op = harmonic_simulation._model.operator("amplitude_fc")
    #     amplitude_op.connect(0, op.outputs.fields_container)
    #     field_ref = amplitude_op.eval()[0]
    #
    #     assert field.component_count == 1
    #     assert field.data.shape == (3,)
    #     assert np.allclose(field.data, field_ref.data)

    # def test_velocity(self, harmonic_simulation):
    #     print(harmonic_simulation)
    #
    #     result = harmonic_simulation.velocity(components=["X"], node_ids=[2, 3, 4])
    #     assert len(result._fc) == 2
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = harmonic_simulation._model.operator("VX")
    #     time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
    #         1, server=harmonic_simulation._model._server
    #     )
    #     op.connect(0, time_scoping)
    #     mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
    #         [2, 3, 4], server=harmonic_simulation._model._server
    #     )
    #     op.connect(1, mesh_scoping)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 1
    #     assert field.data.shape == (3,)
    #     assert np.allclose(field.data, field_ref.data)
    #
    # def test_acceleration(self, harmonic_simulation):
    #     print(harmonic_simulation)
    #
    #     result = harmonic_simulation.acceleration(components=["X"], node_ids=[2, 3, 4])
    #     assert len(result._fc) == 2
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = harmonic_simulation._model.operator("AX")
    #     time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
    #         1, server=harmonic_simulation._model._server
    #     )
    #     op.connect(0, time_scoping)
    #     mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
    #         [2, 3, 4], server=harmonic_simulation._model._server
    #     )
    #     op.connect(1, mesh_scoping)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 1
    #     assert field.data.shape == (3,)
    #     assert np.allclose(field.data, field_ref.data)

    def test_reaction_force(self, allkindofcomplexity):
        harmonic_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.harmonic_mechanical,
        )
        result = harmonic_simulation.reaction_force(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("RF")
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces(self, allkindofcomplexity):
    #     harmonic_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    #     )
    #     result = harmonic_simulation.element_nodal_forces(set_ids=[1])
    #     assert len(result._fc) == 1
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = harmonic_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)
    #
    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces_nodal(self, allkindofcomplexity):
    #     harmonic_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    #     )
    #     result = harmonic_simulation.element_nodal_forces_nodal(set_ids=[1])
    #     assert len(result._fc) == 3
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = harmonic_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     op.connect(9, post.locations.nodal)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)
    #
    # @pytest.mark.skipif(
    #     not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    #     reason="Available starting DPF 5.0",
    # )
    # def test_element_nodal_forces_elemental(self, allkindofcomplexity):
    #     harmonic_simulation = post.load_simulation(
    #         data_sources=allkindofcomplexity,
    #         simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    #     )
    #     result = harmonic_simulation.element_nodal_forces_elemental(set_ids=[1])
    #     assert len(result._fc) == 3
    #     assert result._fc.get_time_scoping().ids == [1]
    #     field = result._fc[0]
    #     op = harmonic_simulation._model.operator("ENF")
    #     op.inputs.bool_rotate_to_global.connect(False)
    #     op.connect(9, post.locations.elemental)
    #     field_ref = op.eval()[0]
    #     assert field.component_count == 3
    #     assert np.allclose(field.data, field_ref.data)

    def test_stress(self, harmonic_simulation):
        # print(harmonic_simulation)
        result = harmonic_simulation.stress(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_elemental(self, harmonic_simulation):
        result = harmonic_simulation.stress_elemental(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_nodal(self, harmonic_simulation):
        result = harmonic_simulation.stress_nodal(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("SX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal(self, harmonic_simulation):
        result = harmonic_simulation.stress_principal(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("S1")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_nodal(self, harmonic_simulation):
        result = harmonic_simulation.stress_principal_nodal(components=2, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("S2")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_elemental(self, harmonic_simulation):
        result = harmonic_simulation.stress_principal_elemental(
            components=3, set_ids=[1]
        )
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("S3")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises(self, harmonic_simulation):
        result = harmonic_simulation.stress_eqv_von_mises(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_elemental(self, harmonic_simulation):
        result = harmonic_simulation.stress_eqv_von_mises_elemental(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_nodal(self, harmonic_simulation):
        result = harmonic_simulation.stress_eqv_von_mises_nodal(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("S_eqv")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elemental_volume(self, harmonic_simulation):
        result = harmonic_simulation.elemental_volume(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("ENG_VOL")
        field_ref = op.eval()[0]
        # print(field_ref)
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain(components=1, set_ids=1)
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_elemental(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_elemental(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_nodal(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_nodal(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPELX")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_principal(components=1, set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        principal_op = harmonic_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_1()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_nodal(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_principal_nodal(
            components=2, set_ids=[1]
        )
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.nodal)
        principal_op = harmonic_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_2()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_elemental(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_principal_elemental(
            components=3, set_ids=[1]
        )
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental)
        principal_op = harmonic_simulation._model.operator(name="invariants_fc")
        principal_op.connect(0, op.outputs.fields_container)
        field_ref = principal_op.outputs.fields_eig_3()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_eqv_von_mises(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = harmonic_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        field_ref = equivalent_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_nodal(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = harmonic_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = harmonic_simulation._model.operator(name="to_nodal_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_eqv_von_mises_elemental(self, harmonic_simulation):
        result = harmonic_simulation.elastic_strain_eqv_von_mises_elemental(set_ids=[1])
        assert len(result._fc) == 2
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("EPEL")
        time_scoping = dpf.time_freq_scoping_factory.scoping_by_set(
            1, server=harmonic_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, post.locations.elemental_nodal)
        equivalent_op = harmonic_simulation._model.operator(name="eqv_fc")
        equivalent_op.connect(0, op.outputs.fields_container)
        average_op = harmonic_simulation._model.operator(name="to_elemental_fc")
        average_op.connect(0, equivalent_op.outputs.fields_container)
        field_ref = average_op.outputs.fields_container()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_disp_external_layer(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.displacement(set_ids=[1], external_layer=True)
        result_all = harmonic_simulation.displacement(set_ids=[1], external_layer=False)
        assert len(result.index.mesh_index) == 4802
        assert len(result.index.mesh_index) == len(result_all.index.mesh_index)
        assert np.allclose(
            result.select(complex=0).max(axis="node_ids").array,
            [2.76941713e-09, 2.76940199e-09, 4.10914311e-10],
        )
        result = harmonic_simulation.displacement(set_ids=[1], external_layer=[1, 2, 3])
        assert len(result.index.mesh_index) == 44
        assert np.allclose(
            result.select(complex=0).max(axis="node_ids").array,
            [-2.50180428e-09, -2.86357660e-10, 8.61977942e-11],
        )

    def test_stress_external_layer(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.stress_elemental(
            all_sets=True, external_layer=True
        )
        assert len(result.index.mesh_index) == 657
        assert len(result.columns.set_ids) == 1
        result = harmonic_simulation.stress_elemental(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        assert len(result.columns.set_ids) == 1
        result = harmonic_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 44
        assert len(result.columns.set_ids) == 1
        assert np.allclose(
            result.select(complex=0).max(axis="node_ids").array, [123327.7207683]
        )

    def test_stress_external_layer2(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], external_layer=True
        )
        assert len(result.index.mesh_index) == 4802
        assert len(result.columns.set_ids) == 1

    def test_strain_external_layer(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.stress_principal_elemental(
            all_sets=True, external_layer=True
        )
        assert len(result.index.mesh_index) == 657
        assert len(result.columns.set_ids) == 1
        assert np.allclose(
            result.select(complex=0).select(set_ids=[1]).max(axis="element_ids").array,
            [1915.12412375],
        )
        result = harmonic_simulation.stress_principal_elemental(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 3
        assert len(result.columns.set_ids) == 1
        result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 44
        assert len(result.columns.set_ids) == 1

    def test_strain_external_layer2(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], external_layer=True
        )
        assert len(result.index.mesh_index) == 4802
        assert len(result.columns.set_ids) == 1

    def test_strain_external_layer3(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.elastic_strain_principal_nodal(
            set_ids=[1], external_layer=[1, 2, 3]
        )
        assert len(result.index.mesh_index) == 44
        assert len(result.columns.set_ids) == 1

    def test_strain_external_layer4(
        self, harmonic_simulation: post.HarmonicMechanicalSimulation
    ):
        result = harmonic_simulation.elastic_strain_eqv_von_mises_elemental(
            set_ids=[1], external_layer=True
        )
        assert len(result.index.mesh_index) == 657
        assert len(result.columns.set_ids) == 1

    def test_disp_skin(self, harmonic_simulation: post.HarmonicMechanicalSimulation):
        if harmonic_simulation._model._server.meet_version("6.2"):
            result = harmonic_simulation.displacement(set_ids=[1], skin=True)
            result_all = harmonic_simulation.displacement(set_ids=[1], skin=False)
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 4184
            else:
                assert len(result.index.mesh_index) == 4802
            assert np.allclose(
                result.select(complex=0).max(axis="node_ids").array,
                [2.76941713e-09, 2.76940199e-09, 4.10914311e-10],
            )
            result = harmonic_simulation.displacement(set_ids=[1], skin=[1, 2, 3])
            assert len(result.index.mesh_index) == 44

    def test_stress_skin(self, harmonic_simulation: post.HarmonicMechanicalSimulation):
        if harmonic_simulation._model._server.meet_version("6.2"):
            result = harmonic_simulation.stress_elemental(all_sets=True, skin=True)
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 1394
            else:
                assert len(result.index.mesh_index) == 3942
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_elemental(
                set_ids=[1], skin=list(range(1, 100))
            )

            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
                assert len(result.index.mesh_index) == 360
            elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 122
            else:
                assert len(result.index.mesh_index) == 192
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_eqv_von_mises_nodal(
                set_ids=[1], skin=list(range(1, 100))
            )

            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
                assert len(result.index.mesh_index) == 1080
            elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 520
            else:
                assert len(result.index.mesh_index) == 530
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_eqv_von_mises_nodal(
                set_ids=[1], skin=True
            )
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 4184
            else:
                assert len(result.index.mesh_index) == 4802
            assert len(result.columns.set_ids) == 1

    def test_strain_skin(self, harmonic_simulation: post.HarmonicMechanicalSimulation):
        if harmonic_simulation._model._server.meet_version("6.2"):
            result = harmonic_simulation.stress_principal_elemental(
                all_sets=True, skin=True
            )
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 1394
            else:
                assert len(result.index.mesh_index) == 3942
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_principal_elemental(
                set_ids=[1], skin=list(range(1, 100))
            )

            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
                assert len(result.index.mesh_index) == 360
            elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 122
            else:
                assert len(result.index.mesh_index) == 192
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(
                set_ids=[1], skin=list(range(1, 100))
            )

            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
                assert len(result.index.mesh_index) == 1080
            elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 520
            else:
                assert len(result.index.mesh_index) == 530
            assert len(result.columns.set_ids) == 1
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
                assert np.allclose(
                    result.select(complex=0).max(axis="node_ids").array,
                    [1.37163319e-06],
                )
            else:
                assert np.allclose(
                    result.select(complex=0).max(axis="node_ids").array,
                    [1.34699501e-06],
                )
            result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(
                set_ids=[1], skin=True
            )
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 4184
            else:
                assert len(result.index.mesh_index) == 4802
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.elastic_strain_principal_nodal(
                set_ids=[1], skin=list(range(1, 100))
            )
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
                assert len(result.index.mesh_index) == 1080
            elif SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 520
            else:
                assert len(result.index.mesh_index) == 530
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.elastic_strain_eqv_von_mises_elemental(
                set_ids=[1],
                skin=True,
            )
            if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1:
                assert len(result.index.mesh_index) == 1394
            else:
                assert len(result.index.mesh_index) == 3942
            assert len(result.columns.set_ids) == 1


def test_elemental_ns_on_nodal_result(modal_frame):
    simulation = post.load_simulation(modal_frame)
    assert "BAR_1" in simulation.named_selections
    disp = simulation.displacement(named_selections=["BAR_1"])
    assert disp.index[0].name == ref_labels.node_ids
    assert len(disp.index[0].values) == 1370


@dataclasses.dataclass
class ReferenceDataItem:
    # Node ids of all the nodes. Note: The list contains duplicates
    # for the nodes at the body interfaces
    node_ids: list[int]
    # The data for each node. The order of the data corresponds to the order of the node ids
    data: list[float]


@dataclasses.dataclass
class ReferenceData:
    # Data Read from Reference CSV files
    # Reference data for all bodies combined
    combined: ReferenceDataItem
    # Reference data per body
    per_id: dict[str, ReferenceDataItem]


def get_node_and_data_map(
    mesh: MeshedRegion, csv_file_path: pathlib.Path
) -> ReferenceDataItem:
    # Get the data from a single csv file and find the corresponding node in the dpf
    # mesh (by coordinates). Note: It looks like the node ids in the csv match the node
    # labels in the dpf mesh. I keep the search for now, but we could probably use the node
    # labels in the csv file directly.
    with open(csv_file_path) as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        next(reader, None)
        node_coordinates = mesh.nodes.coordinates_field
        node_coordinates_csv = []
        data_rows = []
        for idx, row in enumerate(reader):
            node_coordinates_csv.append(
                np.array([float(row[1]), float(row[2]), float(row[3])])
            )
            data_rows.append(float(row[4]))

        node_coordinates_dpf = node_coordinates.data
        node_ids = []
        for row_index, csv_coord in enumerate(node_coordinates_csv):
            index = np.where(
                np.isclose(node_coordinates_dpf, csv_coord, rtol=1e-3).all(axis=1)
            )[0]
            if index.size > 0:
                assert index.size == 1
                node_id = mesh.nodes.scoping.ids[index[0]]

                node_ids.append(node_id)
            else:
                raise RuntimeError(
                    f"Node not found in dpf mesh. Node coordinate: {csv_coord}"
                )

        return ReferenceDataItem(node_ids, data_rows)


def get_ref_data_from_csv(mesh: MeshedRegion, csv_file_name: ReferenceCsvFilesNodal):
    combined_ref_data = get_node_and_data_map(mesh, csv_file_name.combined)
    per_id_ref_data = {}
    for mat_id, csv_file in csv_file_name.per_id.items():
        per_id_ref_data[mat_id] = get_node_and_data_map(mesh, csv_file)
    return ReferenceData(combined_ref_data, per_id_ref_data)


def get_bodies_in_scoping(meshed_region: MeshedRegion, scoping: Scoping):
    elemental_scoping = scoping
    if scoping.location == locations.nodal:
        elemental_scoping = operators.scoping.transpose(
            mesh_scoping=scoping,
            meshed_region=meshed_region,
            inclusive=1,
            requested_location=locations.elemental,
        ).eval()

    mat_field = meshed_region.property_field("mat")
    rescoped_mat_field_op = dpf.operators.scoping.rescope_property_field(
        fields=mat_field,
        mesh_scoping=elemental_scoping,
    )

    rescoped_mat_field = rescoped_mat_field_op.outputs.fields_as_property_field()

    return list(set(rescoped_mat_field.data))


def get_ref_result_per_element(
    csv_file_path: pathlib.Path,
):
    elemental_data = {}

    with open(csv_file_path) as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")

        next(reader, None)
        for idx, row in enumerate(reader):
            element_id = int(row[0])
            assert elemental_data.get(element_id) is None
            elemental_data[element_id] = float(row[1])
    return elemental_data


def get_ref_result_per_node_and_material(
    mesh: MeshedRegion,
    reference_csv_files: ReferenceCsvFilesNodal,
):
    # Get the reference data from the csv files.
    # Returns a dictionary with node_id and mat_id as nested keys.
    # Goes through the nodes and checks which material id contains the node.
    # The combined data is only used for consistency checks: It is checked
    # that the combined data contains the values of each material.
    ref_data = get_ref_data_from_csv(mesh, reference_csv_files)

    node_id_to_row_index_map_combined = {}
    for idx, node_id in enumerate(ref_data.combined.node_ids):
        if node_id is not None:
            if node_id in node_id_to_row_index_map_combined:
                node_id_to_row_index_map_combined[node_id].append(idx)
            else:
                node_id_to_row_index_map_combined[node_id] = [idx]

    # Ensure we have found each node in the input mesh
    assert sorted(node_id_to_row_index_map_combined.keys()) == sorted(
        mesh.nodes.scoping.ids
    )

    data_per_node_and_material = {}

    for node_id, combined_row_indices in node_id_to_row_index_map_combined.items():
        multiplicity_of_node = len(combined_row_indices)
        material_wise_data = {}
        for mat_id, ref_data_item in ref_data.per_id.items():
            if node_id in ref_data_item.node_ids:
                row_index = ref_data_item.node_ids.index(node_id)
                material_wise_data[mat_id] = ref_data_item.data[row_index]

        if len(material_wise_data) != multiplicity_of_node:
            raise RuntimeError(
                f"Inconsistent combined and per material data for node id: {node_id}"
                f" number of entries in combined data: {multiplicity_of_node}, "
                f" number of entries in per material data: {len(material_wise_data)}"
            )

        for mat_id, data_per_material in material_wise_data.items():
            # Check that the data per material is close to one of the
            # data for one of the combined nodes
            assert np.isclose(
                data_per_material,
                np.array(ref_data.combined.data)[combined_row_indices],
            ).any(), f"{node_id}, {mat_id}"

        data_per_node_and_material[node_id] = material_wise_data

    return data_per_node_and_material


def get_ref_per_body_results_mechanical(
    reference_csv_files: ReferenceCsvFilesNodal,
    mesh: MeshedRegion,
):
    return get_ref_result_per_node_and_material(mesh, reference_csv_files)


def get_per_body_results_solid(
    simulation: StaticMechanicalSimulation,
    result_type: str,
    mat_ids: list[int],
    components: list[str],
    additional_scoping: Optional[Scoping],
):
    if additional_scoping:
        transpose_scoping = operators.scoping.transpose()
        transpose_scoping.inputs.mesh_scoping(additional_scoping)
        transpose_scoping.inputs.meshed_region(simulation.mesh._meshed_region)
        transpose_scoping.inputs.inclusive(1)
        transpose_scoping.inputs.requested_location(locations.elemental)

        elemental_scoping = transpose_scoping.eval()

    # Split the mesh by bodies to get an elemental scoping.
    mesh = simulation.mesh._meshed_region
    split_by_property_op = operators.scoping.split_on_property_type()
    split_by_property_op.inputs.mesh(mesh)
    split_by_property_op.inputs.label1("mat")

    body_scopings = split_by_property_op.eval()
    elemental_nodal_result = getattr(simulation, result_type)(components=components)._fc

    solid_mesh = elemental_nodal_result[0].meshed_region
    all_values = {}

    for mat_id in mat_ids:
        body_scoping = body_scopings.get_scoping({"mat": mat_id})
        assert body_scoping.location == locations.elemental

        if additional_scoping is not None:
            scoping_intersect_op = dpf.operators.scoping.intersect()
            scoping_intersect_op.inputs.scopingA.connect(body_scoping)
            scoping_intersect_op.inputs.scopingB.connect(elemental_scoping)

            intersected_scoping = scoping_intersect_op.eval()
            if len(intersected_scoping.ids) == 0:
                continue
        else:
            intersected_scoping = body_scoping

        # Rescope the elemental nodal results to the body
        # and the optional additional scoping
        rescope_op = operators.scoping.rescope_fc()
        rescope_op.inputs.mesh_scoping(intersected_scoping)
        rescope_op.inputs.fields_container(elemental_nodal_result)

        to_nodal_op = dpf.operators.averaging.to_nodal_fc()
        to_nodal_op.inputs.fields_container(rescope_op.outputs.fields_container)

        nodal_fc = None
        if additional_scoping.location == locations.nodal:
            rescope_nodal_op = operators.scoping.rescope_fc()
            rescope_nodal_op.inputs.fields_container(
                to_nodal_op.outputs.fields_container()
            )
            rescope_nodal_op.inputs.mesh_scoping(additional_scoping)
            nodal_fc = rescope_nodal_op.eval()
        else:
            nodal_fc = to_nodal_op.outputs.fields_container()
        assert len(nodal_fc) == 1
        nodal_field = nodal_fc[0]

        values_per_mat = {}
        for node_id in nodal_field.scoping.ids:
            entity_data = nodal_field.get_entity_data_by_id(node_id)
            assert len(entity_data) == 1
            values_per_mat[node_id] = entity_data[0]

        all_values[mat_id] = values_per_mat

    # Get all node_ids so it is easy to build
    # the dictionary with nested labels [node_id][mat_id]
    all_node_ids = set()
    for mat_id in mat_ids:
        all_node_ids.update(all_values[mat_id].keys())

    # Build nested dictionary with node_id and mat_id as nested keys.
    expected_results = {}
    for node_id in all_node_ids:
        expected_results_per_node = {}
        for mat_id in mat_ids:
            if node_id in all_values[mat_id]:
                expected_results_per_node[mat_id] = all_values[mat_id][node_id]
        expected_results[node_id] = expected_results_per_node
    return expected_results


def get_ref_per_body_results_skin(
    simulation: StaticMechanicalSimulation,
    result_type: str,
    mat_ids: list[int],
    components: list[str],
    skin_mesh: MeshedRegion,
    additional_scoping: Optional[Scoping],
):
    if additional_scoping:
        transpose_scoping = operators.scoping.transpose()
        transpose_scoping.inputs.mesh_scoping(additional_scoping)
        transpose_scoping.inputs.meshed_region(simulation.mesh._meshed_region)
        transpose_scoping.inputs.inclusive(1)
        transpose_scoping.inputs.requested_location(locations.elemental)

        elemental_scoping = transpose_scoping.eval()

    # Get the reference skin results.
    # Rescope the skin mesh to each body and extract the corresponding results.

    # Split the mesh by bodies to get an elemental scoping.
    mesh = simulation.mesh._meshed_region
    split_by_property_op = operators.scoping.split_on_property_type()
    split_by_property_op.inputs.mesh(mesh)
    split_by_property_op.inputs.label1("mat")

    body_scopings = split_by_property_op.eval()
    elemental_nodal_result = getattr(simulation, result_type)(components=components)._fc

    skin_values = {}

    solid_mesh = elemental_nodal_result[0].meshed_region

    for mat_id in mat_ids:
        body_scoping = body_scopings.get_scoping({"mat": mat_id})
        assert body_scoping.location == locations.elemental

        if additional_scoping is not None:
            scoping_intersect_op = (
                dpf.operators.scoping.intersect()
            )  # operator instantiation
            scoping_intersect_op.inputs.scopingA.connect(body_scoping)
            scoping_intersect_op.inputs.scopingB.connect(elemental_scoping)

            intersected_scoping = scoping_intersect_op.eval()
            if len(intersected_scoping.ids) == 0:
                continue
        else:
            intersected_scoping = body_scoping

        # Rescope the elemental nodal results to the body
        # The elemental nodal results are used later to get the nodal
        # results
        rescope_op = operators.scoping.rescope_fc()
        rescope_op.inputs.mesh_scoping(intersected_scoping)
        rescope_op.inputs.fields_container(elemental_nodal_result)

        # Rescope the solid mesh
        rescope_mesh_op_solid = operators.mesh.from_scoping()
        rescope_mesh_op_solid.inputs.mesh(solid_mesh)
        rescope_mesh_op_solid.inputs.scoping(intersected_scoping)

        rescoped_solid_mesh = rescope_mesh_op_solid.eval()

        # Get the nodal scoping, which is needed to rescope
        # the skin mesh.
        transpose_scoping = operators.scoping.transpose()
        transpose_scoping.inputs.mesh_scoping(intersected_scoping)
        transpose_scoping.inputs.meshed_region(solid_mesh)
        transpose_scoping.inputs.inclusive(1)

        nodal_scoping = transpose_scoping.eval()

        # Rescope the skin mesh
        rescope_mesh_op_skin = operators.mesh.from_scoping()
        rescope_mesh_op_skin.inputs.mesh(skin_mesh)
        rescope_mesh_op_skin.inputs.scoping(nodal_scoping)
        rescope_mesh_op_skin.inputs.inclusive(0)

        for field in elemental_nodal_result:
            field.meshed_region = rescoped_solid_mesh

        nodal_field = get_expected_nodal_skin_results(
            simulation=simulation,
            result_name=result_type,
            mode=None,
            skin_mesh=rescope_mesh_op_skin.eval(),
            expand_cyclic=False,
            elemental_nodal_results=rescope_op.eval(),
        )

        if additional_scoping and additional_scoping.location == locations.nodal:
            rescope_to_add_scope = operators.scoping.rescope()
            rescope_to_add_scope.inputs.mesh_scoping(additional_scoping)
            rescope_to_add_scope.inputs.fields(nodal_field)
            nodal_field = rescope_to_add_scope.outputs.fields_as_field()

        skin_values_per_mat = {}
        for node_id in nodal_field.scoping.ids:
            entity_data = nodal_field.get_entity_data_by_id(node_id)
            assert len(entity_data) == 1
            skin_values_per_mat[node_id] = entity_data[0]

        skin_values[mat_id] = skin_values_per_mat

    # Get all node_ids so it is easy to build
    # the dictionary with nested labels [node_id][mat_id]
    all_node_ids = set()
    for mat_id in mat_ids:
        all_node_ids.update(skin_values[mat_id].keys())

    # Build nested dictionary with node_id and mat_id as nested keys.
    expected_results = {}
    for node_id in all_node_ids:
        expected_results_per_node = {}
        for mat_id in mat_ids:
            if node_id in skin_values[mat_id]:
                expected_results_per_node[mat_id] = skin_values[mat_id][node_id]
        expected_results[node_id] = expected_results_per_node
    return expected_results


default_per_body_averaging_config = AveragingConfig(
    body_defining_properties=[
        elemental_properties.material,
        "mapdl_element_type_id",
    ],
    average_per_body=True,
)


def get_custom_scope(selection_name: str, mesh: MeshedRegion):
    if selection_name == "BODY_BY_ELEMENT_IDS":
        # Element scope that corresponds to one body
        element_scope = [25, 26, 32, 31, 27, 28, 33, 34, 29, 30, 35, 36]
        custom_scoping = Scoping(ids=element_scope, location=locations.elemental)
        transpose_op = operators.scoping.transpose()
        transpose_op.inputs.requested_location(locations.nodal)
        transpose_op.inputs.inclusive(0)
        transpose_op.inputs.mesh_scoping(custom_scoping)
        transpose_op.inputs.meshed_region(mesh)
        expected_nodal_scope = transpose_op.eval().ids
        return custom_scoping, expected_nodal_scope
    elif selection_name == "SINGLE_NODE":
        expected_nodal_scope = [1]
        custom_scoping = Scoping(ids=expected_nodal_scope, location=locations.nodal)
        return custom_scoping, expected_nodal_scope
    else:
        named_selection_scope = mesh.named_selection("SELECTION")
        assert named_selection_scope.location == locations.nodal
        expected_nodal_scope = named_selection_scope.ids
        transpose_op = operators.scoping.transpose()
        transpose_op.inputs.requested_location(locations.elemental)
        transpose_op.inputs.inclusive(0)
        transpose_op.inputs.mesh_scoping(named_selection_scope)
        transpose_op.inputs.meshed_region(mesh)
        custom_elemental_scoping = transpose_op.eval()

        if selection_name == "SELECTION_CONVERT_TO_ELEMENTAL":
            custom_scoping = Scoping(
                ids=custom_elemental_scoping.ids, location=locations.elemental
            )

        if selection_name == "SELECTION_CONVERT_TO_NODAL":
            custom_scoping = named_selection_scope
        return custom_scoping, expected_nodal_scope


@pytest.mark.parametrize("is_skin", [False, True])
# Note: Selections are only tested on the more complex model (average_per_body_complex_multi_body)
@pytest.mark.parametrize(
    "selection_name",
    [
        None,
        # Use the named selection (nodal selection) in the model to do the selection.
        "SELECTION",
        # Use a custom selection (based on element ids) to do the selection.
        # Selection coincides with one of the bodies
        "BODY_BY_ELEMENT_IDS",
        # Selection of a single node
        "SINGLE_NODE",
        # Use the named selection (nodal selection) in the model, but convert it to
        # node_ids to test the node_ids argument of the results api.
        "SELECTION_CONVERT_TO_NODAL",
        # Use the named selection (nodal selection) in the model, but convert it to
        # element_ids to test the element_ids argument of the results api.
        "SELECTION_CONVERT_TO_ELEMENTAL",
    ],
)
@pytest.mark.parametrize("result", ["stress", "elastic_strain"])
@pytest.mark.parametrize(
    "result_file_str, ref_files",
    [
        (r"average_per_body_two_cubes", "average_per_body_two_cubes_ref"),
        (
            r"average_per_body_complex_multi_body",
            "average_per_body_complex_multi_body_ref",
        ),
    ],
)
def test_averaging_per_body_nodal(
    request, is_skin, result, result_file_str, ref_files, selection_name
):
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0:
        # average per body not supported before 9.0
        return

    ref_files = request.getfixturevalue(ref_files)

    result_file = request.getfixturevalue(result_file_str)

    is_named_selection = selection_name not in [
        "BODY_BY_ELEMENT_IDS",
        "SINGLE_NODE",
        "SELECTION_CONVERT_TO_NODAL",
        "SELECTION_CONVERT_TO_ELEMENTAL",
    ]
    simulation: StaticMechanicalSimulation = post.load_simulation(
        data_sources=result_file,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )
    mesh = simulation.mesh._meshed_region

    expected_nodal_scope = None
    if not is_named_selection:
        custom_scoping, expected_nodal_scope = get_custom_scope(selection_name, mesh)

    components = ["XX"]

    kwargs = {}
    if selection_name is not None:
        if not is_named_selection:
            if result_file_str != "average_per_body_complex_multi_body":
                # Test custom selection only with complex case
                return

            if custom_scoping.location == locations.nodal:
                kwargs["node_ids"] = custom_scoping.ids
            else:
                kwargs["element_ids"] = custom_scoping.ids
        else:
            kwargs["named_selections"] = [selection_name]
    res = simulation._get_result(
        base_name=operator_map[result],
        location=locations.nodal,
        category=ResultCategory.matrix,
        skin=is_skin,
        averaging_config=default_per_body_averaging_config,
        components=components,
        **kwargs,
    )

    named_selection = None
    additional_scoping = None
    if selection_name is None:
        mat_field = mesh.property_field("mat")
        bodies_in_selection = list(set(mat_field.data))

    else:
        if not is_named_selection:
            additional_scoping = custom_scoping
        else:
            additional_scoping = mesh.named_selection(selection_name)
            assert additional_scoping.location == "Nodal"
            named_selection = additional_scoping

        # Get only the bodies that are present in the named selection.
        # Only these bodies are present in the dpf result.
        bodies_in_selection = get_bodies_in_scoping(
            meshed_region=simulation.mesh._meshed_region,
            scoping=additional_scoping,
        )

    if is_skin:
        # Compute reference data on skin (by rescoping results on skin)
        ref_data = get_ref_per_body_results_skin(
            simulation=simulation,
            result_type=result,
            mat_ids=bodies_in_selection,
            components=components,
            skin_mesh=res._fc[0].meshed_region,
            additional_scoping=additional_scoping,
        )
    else:
        # Cannot take reference for Mechanical because the named selection
        # splits a body and therefore the values at the boundaries
        # of the named selection are not the same as in Mechanical
        # Instead the elemental nodal data is rescoped to the additional_scoping and
        # then averaged on that scoping.
        if named_selection is not None or not is_named_selection:
            ref_data = get_per_body_results_solid(
                simulation=simulation,
                result_type=result,
                mat_ids=bodies_in_selection,
                components=components,
                additional_scoping=additional_scoping,
            )
        else:
            # get reference data from mechanical
            ref_data = get_ref_per_body_results_mechanical(ref_files[result], mesh)

    def get_expected_label_space_by_mat_id(mat_id: int):
        # mapdl_element_type_id is not part of the label space before DPF 9.1
        if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1:
            return {
                elemental_properties.material: mat_id,
                "time": 1,
            }
        return {
            elemental_properties.material: mat_id,
            "mapdl_element_type_id": mat_id,
            "time": 1,
        }

    label_spaces_by_mat_id = {}
    for idx in range(len(bodies_in_selection)):
        label_space = res._fc.get_label_space(idx)
        label_spaces_by_mat_id[label_space["mat"]] = label_space

    assert len(label_spaces_by_mat_id) == len(bodies_in_selection)
    for mat_id in bodies_in_selection:
        assert label_spaces_by_mat_id[mat_id] == get_expected_label_space_by_mat_id(
            mat_id
        )

    assert res._fc.get_label_space(len(bodies_in_selection)) == {}

    for node_id in ref_data:
        for mat_id in ref_data[node_id]:
            mat_id_int = int(mat_id)

            if selection_name is not None and mat_id_int not in bodies_in_selection:
                continue
            field = res._fc.get_field({"mat": mat_id_int})

            nodal_value = None
            if expected_nodal_scope is not None:
                assert set(field.scoping.ids).issubset(set(expected_nodal_scope)), set(
                    field.scoping.ids
                ).difference(set(expected_nodal_scope))

                if node_id in expected_nodal_scope:
                    nodal_value = field.get_entity_data_by_id(node_id)
            else:
                nodal_value = field.get_entity_data_by_id(node_id)

            if nodal_value is not None:
                assert np.isclose(
                    nodal_value[0], ref_data[node_id][mat_id], rtol=1e-3
                ), f"{result}, {mat_id}, {node_id}"


@pytest.mark.parametrize("is_skin", [False, True])
@pytest.mark.parametrize("named_selection_name", [None, "SELECTION"])
@pytest.mark.parametrize("result", ["stress", "elastic_strain"])
@pytest.mark.parametrize(
    "result_file",
    [
        r"average_per_body_two_cubes",
        r"average_per_body_complex_multi_body",
    ],
)
def test_averaging_per_body_elemental(
    request, is_skin, result, result_file, named_selection_name
):
    # Expectation is that elemental results are not affected by the average per body flag.

    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0:
        # average per body not supported before 9.0
        return

    result_file = request.getfixturevalue(result_file)
    rst_file = pathlib.Path(result_file)
    simulation: StaticMechanicalSimulation = post.load_simulation(
        data_sources=rst_file,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )

    components = ["XX"]

    named_selections = None
    if named_selection_name is not None:
        named_selections = [named_selection_name]

    kwargs = {
        "base_name": operator_map[result],
        "location": locations.elemental,
        "category": ResultCategory.matrix,
        "skin": is_skin,
        "components": components,
        "named_selections": named_selections,
    }
    res_per_body_fc = simulation._get_result(
        **kwargs, averaging_config=default_per_body_averaging_config
    )._fc

    res_across_bodies_fc = simulation._get_result(
        **kwargs, averaging_config=AveragingConfig()
    )._fc

    assert len(res_across_bodies_fc) == 1
    res_across_bodies_field = res_across_bodies_fc[0]

    mat_property_field = res_across_bodies_field.meshed_region.property_field("mat")
    for element_id in res_across_bodies_field.scoping.ids:
        mat_id_arr = mat_property_field.get_entity_data_by_id(element_id)
        assert len(mat_id_arr) == 1

        res_per_body_field = res_per_body_fc.get_field({"mat": mat_id_arr[0]})
        assert res_across_bodies_field.get_entity_data_by_id(
            element_id
        ) == res_per_body_field.get_entity_data_by_id(element_id)


@pytest.mark.parametrize("is_skin", [False, True])
@pytest.mark.parametrize("average_per_body", [False, True])
@pytest.mark.parametrize("requested_location", ["Nodal", "Elemental"])
def test_build_selection(
    average_per_body_complex_multi_body, average_per_body, is_skin, requested_location
):
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0:
        # Logic has changed with server 9.0
        return

    rst_file = pathlib.Path(average_per_body_complex_multi_body)
    simulation: StaticMechanicalSimulation = post.load_simulation(
        data_sources=rst_file,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )

    scoping = Scoping(
        location=locations.elemental,
        ids=[25, 26, 32, 31, 27, 28, 33, 34, 29, 30, 35, 36],
    )

    selection, rescoping = simulation._build_selection(
        base_name="S",
        category=ResultCategory.matrix,
        location=requested_location,
        skin=is_skin,
        average_per_body=average_per_body,
        selection=None,
        set_ids=None,
        times=None,
        all_sets=True,
        element_ids=scoping.ids,
    )
    selection_wf = selection.spatial_selection._selection
    if selection.spatial_selection.requires_mesh:
        selection_wf.connect(_WfNames.initial_mesh, simulation.mesh._meshed_region)
    scoping_from_selection = selection_wf.get_output(_WfNames.scoping, Scoping)

    if is_skin or average_per_body:
        # If request is for skin or average per body, the location should be elemental
        # because force_elemental_nodal is True
        assert scoping_from_selection.location == locations.elemental
        assert set(scoping_from_selection.ids) == set(scoping.ids)
    else:
        assert scoping_from_selection.location == requested_location
        if requested_location == locations.nodal:
            assert len(scoping_from_selection.ids) == 36
        else:
            assert set(scoping_from_selection.ids) == set(scoping.ids)


def test_beam_results_on_skin(beam_example):
    simulation: StaticMechanicalSimulation = post.load_simulation(
        data_sources=beam_example,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_10_0:
        # Add beams on skin not activated before 10.0
        return

    res = simulation.displacement(skin=True, norm=True)

    element_type_array = res._fc[0].meshed_region.elements.element_types_field.data
    element_count_dict = {
        key: sum(1 for _ in value) for key, value in groupby(sorted(element_type_array))
    }

    unit_converter = dpf.operators.math.unit_convert(
        unit_name=2,  # NMM unit system
    )

    unit_converter.inputs.entity_to_convert(res._fc[0])
    converted_field = unit_converter.eval()

    assert element_types.Line2.value in element_count_dict.keys()

    assert element_count_dict[element_types.Line2.value] == 40

    assert converted_field.max().data[0] == pytest.approx(190, 1e-2)
