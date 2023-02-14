import ansys.dpf.core as core
import numpy as np
import pytest
from pytest import fixture

import ansys.dpf.post as dpf
from ansys.dpf.post.common import AvailableSimulationTypes


@fixture
def static_simulation(static_rst):
    return dpf.load_simulation(
        data_sources=static_rst,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )


@fixture
def transient_simulation(plate_msup):
    return dpf.load_simulation(
        data_sources=plate_msup,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )


def test_simulation_results(static_simulation):
    results = static_simulation.results
    assert len(results) == 12
    assert all(isinstance(x, str) for x in results)


def test_simulation_geometries(static_simulation):
    geometries = static_simulation.geometries
    assert geometries == []


def test_simulation_boundary_conditions(static_simulation):
    boundary_conditions = static_simulation.boundary_conditions
    assert boundary_conditions == []


def test_simulation_loads(static_simulation):
    loads = static_simulation.loads
    assert loads == []


def test_simulation_mesh(static_simulation):
    mesh = static_simulation.mesh
    assert isinstance(mesh, dpf.mesh.Mesh)


def test_simulation_named_selections(static_simulation):
    named_selections = static_simulation.named_selections
    assert len(named_selections) == 1
    assert all(isinstance(x, str) for x in named_selections)


def test_simulation_active_selection(static_simulation):
    assert static_simulation.active_selection is None
    selection = dpf.selection.Selection()
    static_simulation.activate_selection(selection=selection)
    assert static_simulation.active_selection == selection
    static_simulation.deactivate_selection()
    assert static_simulation.active_selection is None


def test_simulation_plot(static_simulation):
    static_simulation.plot()


class TestStaticMechanicalSimulation:
    def test_warning_empty(self, static_simulation):
        with pytest.warns(expected_warning=UserWarning, match="empty"):
            _ = static_simulation.displacement(
                components=1, node_ids=[1001, 1002, 1003]
            )

    def test_displacement(self, static_simulation):
        displacement_x = static_simulation.displacement(
            components=["X"], node_ids=[42, 43, 44], set_ids=[1]
        )
        assert len(displacement_x._fc) == 1
        assert displacement_x._fc.get_time_scoping().ids == [1]
        field = displacement_x._fc[0]
        op = static_simulation._model.operator("UX")
        mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
            [42, 43, 44], server=static_simulation._model._server
        )
        op.connect(1, mesh_scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(field.data, field_ref.data)

        displacement_y = static_simulation.displacement(
            components=["2"],
            named_selections=static_simulation.named_selections[0],
            load_steps=(1, 1),
        )
        assert len(displacement_y._fc) == 1
        assert displacement_y._fc.get_time_scoping().ids == [1]
        field = displacement_y._fc[0]
        op = static_simulation._model.operator("UY")
        mesh_scoping = core.mesh_scoping_factory.named_selection_scoping(
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
        mesh_scoping = core.mesh_scoping_factory.named_selection_scoping(
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
        mesh_scoping = core.mesh_scoping_factory.elemental_scoping(
            element_ids=[1, 2, 3],
            server=static_simulation._model._server,
        )
        mesh_scoping = core.operators.scoping.transpose(
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
        mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
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
        op.connect(9, core.locations.elemental_nodal)
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
        op.connect(9, core.locations.elemental)
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
        op.connect(9, core.locations.nodal)
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
        op.connect(9, core.locations.elemental_nodal)
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
        op.connect(9, core.locations.nodal)
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
        op.connect(9, core.locations.elemental)
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
        op.connect(9, core.locations.elemental_nodal)
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
        op.connect(9, core.locations.elemental)
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
        op.connect(9, core.locations.nodal)
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
        op.connect(9, core.locations.nodal)
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
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (12,)
        assert np.allclose(field.data, field_ref.data)

    def test_element_nodal_forces(self, allkindofcomplexity):
        static_simulation = dpf.load_simulation(data_sources=allkindofcomplexity)
        element_nodal_forces = static_simulation.element_nodal_forces()
        assert len(element_nodal_forces._fc) == 1
        assert element_nodal_forces._fc.get_time_scoping().ids == [1]
        field = element_nodal_forces._fc[0]
        op = static_simulation._model.operator("ENF")
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert field.data.shape == (103766, 3)
        assert np.allclose(field.data, field_ref.data)

    def test_element_nodal_forces_nodal(self, allkindofcomplexity):
        static_simulation = dpf.load_simulation(data_sources=allkindofcomplexity)
        element_nodal_forces = static_simulation.element_nodal_forces_nodal()
        assert len(element_nodal_forces._fc) == 3
        assert element_nodal_forces._fc.get_time_scoping().ids == [1]
        field = element_nodal_forces._fc[0]
        op = static_simulation._model.operator("ENF")
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert field.data.shape == (14982, 3)
        assert np.allclose(field.data, field_ref.data)

    def test_element_nodal_forces_elemental(self, allkindofcomplexity):
        static_simulation = dpf.load_simulation(data_sources=allkindofcomplexity)
        element_nodal_forces = static_simulation.element_nodal_forces_elemental()
        assert len(element_nodal_forces._fc) == 3
        assert element_nodal_forces._fc.get_time_scoping().ids == [1]
        field = element_nodal_forces._fc[0]
        op = static_simulation._model.operator("ENF")
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert field.data.shape == (9433, 3)
        assert np.allclose(field.data, field_ref.data)


class TestTransientMechanicalSimulation:
    def test_displacement(self, transient_simulation):
        result = transient_simulation.displacement(
            components=["X"], node_ids=[2, 3, 4], set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("UX")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
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
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
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
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
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
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_elemental(self, transient_simulation):
        result = transient_simulation.stress_elemental(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("SX")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_nodal(self, transient_simulation):
        result = transient_simulation.stress_nodal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("SX")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal(self, transient_simulation):
        result = transient_simulation.stress_principal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S1")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_principal_nodal(self, transient_simulation):
        result = transient_simulation.stress_principal_nodal(components=2, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S2")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.nodal)
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
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises(self, transient_simulation):
        result = transient_simulation.stress_eqv_von_mises(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S_eqv")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_elemental(self, transient_simulation):
        result = transient_simulation.stress_eqv_von_mises_elemental(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S_eqv")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_nodal(self, transient_simulation):
        result = transient_simulation.stress_eqv_von_mises_nodal(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("S_eqv")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain(self, transient_simulation):
        result = transient_simulation.elastic_strain(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPELX")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental_nodal)
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
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_nodal(self, transient_simulation):
        result = transient_simulation.elastic_strain_nodal(components=1, set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPELX")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.nodal)
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
        op = transient_simulation._model.operator("EPEL1")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental_nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_nodal(self, transient_simulation):
        result = transient_simulation.elastic_strain_principal_nodal(
            components=2, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL2")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_elastic_strain_principal_elemental(self, transient_simulation):
        result = transient_simulation.elastic_strain_principal_elemental(
            components=3, set_ids=[2]
        )
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("EPEL3")
        time_scoping = core.time_freq_scoping_factory.scoping_by_set(
            2, server=transient_simulation._model._server
        )
        op.connect(0, time_scoping)
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_reaction_force(self, allkindofcomplexity):
        transient_simulation = dpf.load_simulation(
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
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_structural_temperature_elemental(self, transient_simulation):
        result = transient_simulation.structural_temperature_elemental(set_ids=[2])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [2]
        field = result._fc[0]
        op = transient_simulation._model.operator("BFE")
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert np.allclose(field.data, field_ref.data)

    def test_element_nodal_forces(self, allkindofcomplexity):
        transient_simulation = dpf.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.element_nodal_forces(set_ids=[1])
        assert len(result._fc) == 1
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENF")
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    def test_element_nodal_forces_nodal(self, allkindofcomplexity):
        transient_simulation = dpf.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.element_nodal_forces_nodal(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENF")
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)

    def test_element_nodal_forces_elemental(self, allkindofcomplexity):
        transient_simulation = dpf.load_simulation(
            data_sources=allkindofcomplexity,
            simulation_type=AvailableSimulationTypes.transient_mechanical,
        )
        result = transient_simulation.element_nodal_forces_elemental(set_ids=[1])
        assert len(result._fc) == 3
        assert result._fc.get_time_scoping().ids == [1]
        field = result._fc[0]
        op = transient_simulation._model.operator("ENF")
        op.connect(9, core.locations.elemental)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        assert np.allclose(field.data, field_ref.data)
