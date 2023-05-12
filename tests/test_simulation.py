import os.path

import ansys.dpf.core as dpf
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_2,
)
import numpy as np
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post.common import AvailableSimulationTypes
from ansys.dpf.post.index import ref_labels


@fixture
def static_simulation(static_rst):
    return post.load_simulation(
        data_sources=static_rst,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )


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


def test_simulation_results(static_simulation):
    results = static_simulation.results
    assert len(results) == 12
    assert all(isinstance(x, str) for x in results)


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
    static_simulation.plot()


class TestStaticMechanicalSimulation:
    def test_cyclic(self, simple_cyclic):
        simulation = post.StaticMechanicalSimulation(simple_cyclic)
        result = simulation.stress(expand_cyclic=False)
        print(result)
        assert "base_sector" in result.columns.names
        result = simulation.stress(expand_cyclic=True)
        print(result)
        assert "base_sector" not in result.columns.names

    def test_multi_stage(self, multi_stage_cyclic):
        simulation = post.StaticMechanicalSimulation(multi_stage_cyclic)
        result = simulation.stress(expand_cyclic=False)
        print(result)
        assert "base_sector" in result.columns.names
        assert "stage" in result.columns.names
        result = simulation.stress(expand_cyclic=True)
        print(result)
        assert "base_sector" not in result.columns.names
        assert "stage" not in result.columns.names

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
        print(field_ref)
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

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces(self, allkindofcomplexity):
        static_simulation = post.load_simulation(data_sources=allkindofcomplexity)
        element_nodal_forces = static_simulation.element_nodal_forces()
        assert len(element_nodal_forces._fc) == 1
        assert element_nodal_forces._fc.get_time_scoping().ids == [1]
        field = element_nodal_forces._fc[0]
        op = static_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert field.data.shape == (103750, 3)
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_nodal(self, allkindofcomplexity):
        static_simulation = post.load_simulation(data_sources=allkindofcomplexity)
        element_nodal_forces = static_simulation.element_nodal_forces_nodal()
        assert len(element_nodal_forces._fc) == 3
        assert element_nodal_forces._fc.get_time_scoping().ids == [1]
        field = element_nodal_forces._fc[0]
        op = static_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert field.data.shape == (14982, 3)
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_elemental(self, allkindofcomplexity):
        static_simulation = post.load_simulation(data_sources=allkindofcomplexity)
        element_nodal_forces = static_simulation.element_nodal_forces_elemental()
        assert len(element_nodal_forces._fc) == 3
        assert element_nodal_forces._fc.get_time_scoping().ids == [1]
        field = element_nodal_forces._fc[0]
        op = static_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert field.data.shape == (9433, 3)
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
        assert len(result.index.mesh_index) == 18

    def test_skin_layer4(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.stress_principal_nodal(skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 44

    def test_skin_layer5(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.elastic_strain_eqv_von_mises_nodal(skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 44

    def test_skin_layer6(self, static_simulation: post.StaticMechanicalSimulation):
        result = static_simulation.stress_principal_elemental(skin=[1, 2, 3])
        assert len(result.index.mesh_index) == 18


class TestTransientMechanicalSimulation:
    @fixture
    def transient_simulation(self, plate_msup):
        return post.load_simulation(
            data_sources=plate_msup,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )

    def test_times_argument(self, transient_simulation, static_simulation):
        with pytest.raises(
            ValueError, match="Could not find time=0.0s in the simulation."
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
            ValueError, match="Could not find time=0.1495s in the simulation."
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
        print(field_ref)
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

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces(self, allkindofcomplexity):
        transient_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.element_nodal_forces(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_nodal(self, allkindofcomplexity):
        transient_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.element_nodal_forces_nodal(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_elemental(self, allkindofcomplexity):
        transient_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.element_nodal_forces_elemental(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

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
        assert len(result.index.mesh_index) == 240

    def test_skin_layer4(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.elastic_strain_eqv_von_mises_elemental(
            skin=list(range(1, 100))
        )
        assert len(result.index.mesh_index) == 240

    def test_skin_layer5(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.stress_principal_nodal(skin=list(range(1, 100)))
        assert len(result.index.mesh_index) == 393

    def test_skin_layer6(
        self, transient_simulation: post.TransientMechanicalSimulation
    ):
        result = transient_simulation.elastic_strain_eqv_von_mises_nodal(
            skin=list(range(1, 100))
        )
        assert len(result.index.mesh_index) == 393


class TestModalMechanicalSimulation:
    @fixture
    def modal_simulation(self, modalallkindofcomplexity):
        return post.load_simulation(
            data_sources=modalallkindofcomplexity,
            simulation_type=AvailableSimulationTypes.modal_mechanical,
        )

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
        print(displacement)
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
        print(modal_simulation)
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

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces(self, allkindofcomplexity):
        modal_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.modal_mechanical,
        )
        result = modal_simulation.element_nodal_forces(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_nodal(self, allkindofcomplexity):
        modal_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.modal_mechanical,
        )
        result = modal_simulation.element_nodal_forces_nodal(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_elemental(self, allkindofcomplexity):
        modal_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.modal_mechanical,
        )
        result = modal_simulation.element_nodal_forces_elemental(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = modal_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

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
        print(field_ref)
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
        if frame_modal_simulation._model._server.meet_version("6.2"):
            result = frame_modal_simulation.stress_elemental(all_sets=True, skin=True)
            assert len(result.index.mesh_index) == 11146
            assert len(result.columns.set_ids) == 6
        result = frame_modal_simulation.stress_elemental(
            set_ids=[1], skin=list(range(1, 100))
        )
        assert len(result.index.mesh_index) == 110
        assert len(result.columns.set_ids) == 1
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
        assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [1285.17926915])
        result = frame_modal_simulation.stress_eqv_von_mises_nodal(
            set_ids=[1], skin=True
        )
        assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1

    def test_strain_skin(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        if frame_modal_simulation._model._server.meet_version("6.2"):
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
        assert len(result.index.mesh_index) == 110
        assert len(result.columns.set_ids) == 1

    def test_strain_skin2(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], skin=frame_modal_simulation.mesh.element_ids
        )
        assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1
        assert np.allclose(result.max(axis="node_ids").array, [0.00684776])

    def test_strain_skin3(self, frame_modal_simulation: post.ModalMechanicalSimulation):
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_nodal(
            set_ids=[1], skin=True
        )
        assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1
        result = frame_modal_simulation.elastic_strain_principal_nodal(
            set_ids=[1], skin=frame_modal_simulation.mesh.element_ids
        )
        assert len(result.index.mesh_index) == 5828
        assert len(result.columns.set_ids) == 1
        result = frame_modal_simulation.elastic_strain_eqv_von_mises_elemental(
            set_ids=[1], skin=True
        )
        assert len(result.index.mesh_index) == 11146
        assert len(result.columns.set_ids) == 1


class TestHarmonicMechanicalSimulation:
    @fixture
    def harmonic_simulation(self, complex_model):
        return post.load_simulation(
            data_sources=complex_model,
            simulation_type=AvailableSimulationTypes.harmonic_mechanical,
        )

    def test_cyclic(self, simple_cyclic):
        simulation = post.HarmonicMechanicalSimulation(simple_cyclic)
        result = simulation.displacement(expand_cyclic=False)
        print(result)
        assert "base_sector" in result.columns.names
        result = simulation.displacement(expand_cyclic=True)
        print(result)
        assert "base_sector" not in result.columns.names

    def test_multi_stage(self, multi_stage_cyclic):
        simulation = post.HarmonicMechanicalSimulation(multi_stage_cyclic)
        result = simulation.displacement(expand_cyclic=False)
        print(result)
        assert "base_sector" in result.columns.names
        assert "stage" in result.columns.names
        result = simulation.displacement(expand_cyclic=True)
        print(result)
        assert "base_sector" not in result.columns.names
        assert "stage" not in result.columns.names

    def test_displacement(self, harmonic_simulation):
        print(harmonic_simulation)

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

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces(self, allkindofcomplexity):
        harmonic_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.harmonic_mechanical,
        )
        result = harmonic_simulation.element_nodal_forces(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_nodal(self, allkindofcomplexity):
        harmonic_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.harmonic_mechanical,
        )
        result = harmonic_simulation.element_nodal_forces_nodal(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    @pytest.mark.skipif(
        not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0,
        reason="Available starting DPF 5.0",
    )
    def test_element_nodal_forces_elemental(self, allkindofcomplexity):
        harmonic_simulation = post.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.harmonic_mechanical,
        )
        result = harmonic_simulation.element_nodal_forces_elemental(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = harmonic_simulation._model.operator("ENF")
        op.inputs.bool_rotate_to_global.connect(False)
        op.connect(9, post.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    def test_stress(self, harmonic_simulation):
        print(harmonic_simulation)
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
        print(field_ref)
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
            assert len(result.index.mesh_index) == 3942
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_elemental(
                set_ids=[1], skin=list(range(1, 100))
            )
            assert len(result.index.mesh_index) == 192
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_eqv_von_mises_nodal(
                set_ids=[1], skin=list(range(1, 100))
            )
            assert len(result.index.mesh_index) == 530
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_eqv_von_mises_nodal(
                set_ids=[1], skin=True
            )
            assert len(result.index.mesh_index) == 4802
            assert len(result.columns.set_ids) == 1

    def test_strain_skin(self, harmonic_simulation: post.HarmonicMechanicalSimulation):
        if harmonic_simulation._model._server.meet_version("6.2"):
            result = harmonic_simulation.stress_principal_elemental(
                all_sets=True, skin=True
            )
            assert len(result.index.mesh_index) == 3942
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.stress_principal_elemental(
                set_ids=[1], skin=list(range(1, 100))
            )
            assert len(result.index.mesh_index) == 192
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(
                set_ids=[1], skin=list(range(1, 100))
            )
            assert len(result.index.mesh_index) == 530
            assert len(result.columns.set_ids) == 1
            assert np.allclose(
                result.select(complex=0).max(axis="node_ids").array, [1.34699501e-06]
            )
            result = harmonic_simulation.elastic_strain_eqv_von_mises_nodal(
                set_ids=[1], skin=True
            )
            assert len(result.index.mesh_index) == 4802
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.elastic_strain_principal_nodal(
                set_ids=[1], skin=list(range(1, 100))
            )
            assert len(result.index.mesh_index) == 530
            assert len(result.columns.set_ids) == 1
            result = harmonic_simulation.elastic_strain_eqv_von_mises_elemental(
                set_ids=[1], skin=True
            )
            assert len(result.index.mesh_index) == 3942
            assert len(result.columns.set_ids) == 1


def test_elemental_ns_on_nodal_result(modal_frame):
    simulation = post.load_simulation(modal_frame)
    assert "BAR_1" in simulation.named_selections
    disp = simulation.displacement(named_selections=["BAR_1"])
    assert disp.index[0].name == ref_labels.node_ids
    assert len(disp.index[0].values) == 1370
