import ansys.dpf.core as core
import numpy as np
from pytest import fixture

import ansys.dpf.post as dpf


@fixture
def static_simulation(static_rst):
    return dpf.load_simulation(data_sources=static_rst)


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
    def test_displacement(self, static_simulation):
        displacement_x = static_simulation.displacement(
            components=["X"], nodes=[42, 43, 44], set_ids=[1]
        )
        assert len(displacement_x._fc) == 1
        assert displacement_x._fc.time_freq_support.time_frequencies.data == 1
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
            named_selection=static_simulation.named_selections[0],
            load_steps=[1],
            sub_steps=[1],
        )
        assert len(displacement_y._fc) == 1
        assert displacement_y._fc.time_freq_support.time_frequencies.data == 1
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
            named_selection=static_simulation.named_selections[0],
            load_steps=1,
            sub_steps=1,
        )
        assert len(displacement_z._fc) == 1
        assert displacement_z._fc.time_freq_support.time_frequencies.data == 1
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
            elements=[1, 2, 3],
            set_ids=1,
        )
        assert len(displacement_z._fc) == 1
        assert displacement_z._fc.time_freq_support.time_frequencies.data == 1
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

    def test_stress(self, static_simulation):
        stress_x = static_simulation.stress(components=1)
        assert len(stress_x._fc) == 1
        assert stress_x._fc.time_freq_support.time_frequencies.data == 1
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
        assert stress_x._fc.time_freq_support.time_frequencies.data == 1
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
        assert stress_x._fc.time_freq_support.time_frequencies.data == 1
        field = stress_x._fc[0]
        op = static_simulation._model.operator("SX")
        op.connect(9, core.locations.nodal)
        field_ref = op.eval()[0]
        assert field.component_count == 1
        assert field.data.shape == (81,)
        assert np.allclose(field.data, field_ref.data)

    def test_stress_eqv_von_mises_elemental(self, static_simulation):
        stress_vm = static_simulation.stress_eqv_von_mises_elemental()
        assert len(stress_vm._fc) == 1
        assert stress_vm._fc.time_freq_support.time_frequencies.data == 1
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
        assert stress_vm._fc.time_freq_support.time_frequencies.data == 1
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
        assert reaction_force._fc.time_freq_support.time_frequencies.data == 1
        field = reaction_force._fc[0]
        op = static_simulation._model.operator("RF")
        op.connect(1, static_simulation.mesh._meshed_region.nodes.scoping)
        field_ref = op.eval()[0]
        assert field.component_count == 3
        # assert field.data.shape == (81, 3)
        assert np.allclose(field.data, field_ref.data)

    # def test_stress_elemental_volume(self, static_simulation):
    #     elemental_volume = static_simulation.elemental_volume()
    #     assert len(elemental_volume._fc) == 1
    #     assert elemental_volume._fc.time_freq_support.time_frequencies.data == 1
    #     field = elemental_volume._fc[0]
    #     op = static_simulation._model.operator("ENG_VOL")
    #     # op.connect(9, core.locations.elemental)
    #     field_ref = op.eval()[0]
    #     print(field_ref)
    #     assert field.component_count == 1
    #     assert field.data.shape == (8,)
    #     assert np.allclose(field.data, field_ref.data)
