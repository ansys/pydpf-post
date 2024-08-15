import os.path
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
)
import numpy as np
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post import StaticMechanicalSimulation
from ansys.dpf.post.common import AvailableSimulationTypes, elemental_properties
from ansys.dpf.post.index import ref_labels
from ansys.dpf.post.meshes import Meshes
from ansys.dpf.post.simulation import Simulation
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_2,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0,
)


def is_principal(mode: str) -> bool:
    return mode == "principal"


def is_equivalent(mode: str) -> bool:
    return mode == "equivalent"


def mode_suffix(mode: str) -> str:
    if mode == "equivalent":
        return "_eqv_von_mises"
    elif mode == "principal":
        return "_principal"
    return ""


def get_expected_average_skin_value(
    element_id: int,
    solid_mesh: MeshedRegion,
    skin_mesh: MeshedRegion,
    elemental_nodal_solid_data_field: Field,
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
        expected_average_skin_values[skin_element_id] = average
    return expected_average_skin_values


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

    for element_id in element_ids:
        expected_skin_values = get_expected_average_skin_value(
            element_id=element_id,
            solid_mesh=fc_elemental_nodal[0].meshed_region,
            skin_mesh=result_skin_scoped_elemental._fc[0].meshed_region,
            elemental_nodal_solid_data_field=fc_elemental_nodal[0],
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
            ) in expected_skin_values.items():
                field.append(list(expected_skin_value), skin_element_id)
            if is_principal(mode):
                invariant_op = operators.invariant.principal_invariants()
                invariant_op.inputs.field(field)
                field_out = invariant_op.outputs.field_eig_1()
            else:
                invariant_op = static_simulation._model.operator(name="eqv")
                invariant_op.inputs.field(field)
                field_out = invariant_op.outputs.field()

            for skin_element_id in expected_skin_values:
                expected_skin_values[skin_element_id] = field_out.get_entity_data_by_id(
                    skin_element_id
                )

        for skin_element_id, expected_skin_value in expected_skin_values.items():
            actual_skin_value = result_skin_scoped_elemental._fc[
                0
            ].get_entity_data_by_id(skin_element_id)
            assert np.allclose(actual_skin_value, expected_skin_value)


def copy_fields_container(fc):
    """
    Create a copy of a fields container. Just works with a single field
    and a symmetric matrix field.
    """
    fields_container = FieldsContainer()
    field = Field(nature=natures.symmatrix)
    for entity_id in fc[0].scoping.ids:
        entity_data = fc[0].get_entity_data_by_id(entity_id)
        field.append(entity_data, entity_id)
    fields_container.add_label("time")
    fields_container.add_field({"time": 1}, field)
    return fields_container


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


def get_nodal_results(
    simulation: Simulation,
    result_name: str,
    mode: str,
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
    # which corresponds to the case 2 above. We should probably fix this, so it
    # matches the data of case 1

    # If we don't get a elemental_nodal_result_op, this means the result does not support
    # elemental nodal evaluation. Currently used only for displacement.
    # In this case we just get the nodal results directly (case 1 above)
    if elemental_nodal_results is None:
        nodal_result_op = simulation._model.operator(name=operator_map[result_name])
        if hasattr(nodal_result_op, "requested_location"):
            # Displacement operator does not have the requested location input
            nodal_result_op.inputs.requested_location("Nodal")
        nodal_result_op.inputs.time_scoping([1])
        fields_container = nodal_result_op.outputs.fields_container()

    else:
        if is_equivalent(mode) and result_name == "elastic_strain":
            # For elastic strain results, the computation of the equivalent
            # value happens before the averaging.
            invariant_op = simulation._model.operator(name="eqv_fc")
            invariant_op.inputs.fields_container(elemental_nodal_results)
            fields_container = invariant_op.outputs.fields_container()
        else:
            fields_container = elemental_nodal_results

    to_nodal = operators.averaging.to_nodal_fc()
    to_nodal.inputs.fields_container(fields_container)
    fc_nodal = to_nodal.outputs.fields_container()

    if is_principal(mode):
        if result_name == "elastic_strain":
            # For the elastic strain result, the invariants_fc operator
            # multiplies the off-diagonal components by 0.5
            # It checks if the field has a "strain" property. See InvariantOperators.cpp
            # line 1024. Since the field properties are not exposed in python
            # we copy the strain field in a new field that does not specify the property.
            # This is just to make the tests pass. The actual solution is
            # to preserve the strain property in the skin workflow.
            # See also test_strain_scaling above.
            fields_container = copy_fields_container(fc_nodal)
        else:
            fields_container = fc_nodal
        invariant_op = simulation._model.operator(name="invariants_fc")
        invariant_op.inputs.fields_container(fields_container)
        fc_nodal = invariant_op.outputs.fields_eig_1()

    if is_equivalent(mode) and result_name != "elastic_strain":
        # The copy is needed for the same reason as above for the
        # principal results.
        # For elastic strain results, the computation of the equivalent
        # value happens before the averaging.
        fields_container = copy_fields_container(fc_nodal)

        invariant_op = simulation._model.operator(name="eqv_fc")
        invariant_op.inputs.fields_container(fields_container)
        fc_nodal = invariant_op.outputs.fields_container()
    return fc_nodal


@fixture
def static_simulation(static_rst):
    return post.load_simulation(
        data_sources=static_rst,
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

    def test_raise_node_ids_elemental(self, static_simulation):
        with pytest.raises(
            ValueError, match="Argument 'node_ids' can only be used if 'location'"
        ):
            _ = static_simulation.stress(
                node_ids=[42], location=post.locations.elemental
            )

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

    def test_strain_scaling(self, static_simulation: post.StaticMechanicalSimulation):
        # Documents that strain scaling is needed to get the same result, when
        # the strain is copied to a new field.
        strain_result = static_simulation.elastic_strain_elemental()

        invariant_op = static_simulation._model.operator(name="invariants_fc")
        invariant_op.inputs.fields_container(strain_result._fc)
        invariant_fc = invariant_op.outputs.fields_eig_1()

        fields_container = FieldsContainer()
        field = Field(nature=natures.symmatrix)
        for entity_id in strain_result._fc[0].scoping.ids:
            entity_data = strain_result._fc[0].get_entity_data_by_id(entity_id)
            # scaling is needed
            entity_data[:, 3:7] = entity_data[:, 3:7] / 2
            field.append(entity_data, entity_id)

        fields_container.add_label("time")
        fields_container.add_field({"time": 1}, field)
        invariant_op.inputs.fields_container(fields_container)
        invariant_fc_with_scale = invariant_op.outputs.fields_eig_1()

        assert invariant_fc[0].get_entity_data_by_id(1) == invariant_fc_with_scale[
            0
        ].get_entity_data_by_id(1)


# List of element configurations for each simulation type
# The commented configurations contain "corners" which yield incorrect
# results because the skin data is interpolated
element_configurations = {
    "static_simulation": {
        1: [1],
        2: [1, 2],
        # 3: [1, 2, 3],
        4: [1, 2, 3, 4],
        # 5: [1, 2, 3, 4, 5],
        # 6: [1, 2, 3, 4, 5, 6],
        # 7: [1, 2, 3, 4, 5, 6, 7],
        8: [1, 2, 3, 4, 5, 6, 7, 8],
    },
    "transient_simulation": {
        1: [1],
        2: [1, 2],
        #  3: [1, 2, 3],
        4: [1, 2, 3, 4],
        # 5: [1, 2, 3, 4, 5],
        # 6: [1, 2, 3, 4, 5, 6],
        # 7: [1, 2, 3, 4, 5, 6, 7],
        8: [1, 2, 3, 4, 5, 6, 7, 8],
    },
    "modal_simulation": {
        1: [1],
        2: [1, 2],
        3: [1, 2, 3],
        #  4: [1, 2, 3, 4, 5, 6, 7, 8]
    },
    "harmonic_simulation": {
        1: [1],
        2: [1, 2],
        3: [1, 2, 3],
        #  4: list(range(1, 100))
    },
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
    if not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
        # Before 8.0, the solid mesh cannot be connected to the solid_to_skin
        # operator. This yield incorrect results. Therefore we skip all the tests
        # for older versions.
        return

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
            solid_elements = simulation.split_mesh_by_properties(
                {elemental_properties.element_type: element_types.Hex20.value}
            )
            element_ids = solid_elements.element_ids
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

    # For displacements the nodal result
    # is just called displacement without
    # the "nodal" suffix
    nodal_suffix = "_nodal"
    if result_name == "displacement":
        nodal_suffix = ""
    if not is_cyclic_simulation:
        fc_nodal = get_nodal_results(
            simulation=simulation,
            result_name=result_name,
            mode=mode,
            elemental_nodal_results=fc_elemental_nodal,
        )
    else:
        fc_nodal = getattr(
            simulation, f"{result_name}{mode_suffix(mode)}{nodal_suffix}"
        )(set_ids=[1], skin=skin, expand_cyclic=True)._fc

    # Not all the simulation types have the expand_cyclic argument
    kwargs = {}
    if is_cyclic_simulation:
        kwargs["expand_cyclic"] = True

    result_skin_scoped_nodal = getattr(
        simulation, f"{result_name}{mode_suffix(mode)}{nodal_suffix}"
    )(set_ids=[1], skin=skin, **kwargs)
    nodal_skin_field = result_skin_scoped_nodal._fc[0]

    for node_id in nodal_skin_field.scoping.ids:
        assert np.allclose(
            fc_nodal[0].get_entity_data_by_id(node_id),
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
