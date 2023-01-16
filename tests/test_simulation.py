import numpy as np
from pytest import fixture

import ansys.dpf.post as dpf


@fixture
def static_simulation(simple_bar):
    return dpf.load_simulation(data_sources=simple_bar)


def test_simulation_results(static_simulation):
    results = static_simulation.results
    assert len(results) == 10
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
            components=["X"], nodes=[1, 2, 3], set_ids=[1]
        )
        assert len(displacement_x._fc) == 1
        assert displacement_x._fc.time_freq_support.time_frequencies.data == 1
        field = displacement_x._fc[0]
        assert field.component_count == 1
        assert field.data.shape == (3,)
        assert np.allclose(
            field.data, [-6.57505989e-08, -3.54701562e-08, -1.44993319e-08]
        )
        displacement_y = static_simulation.displacement(
            components=["2"],
            named_selection=static_simulation.named_selections[0],
            load_steps=[1],
            sub_steps=[1],
        )
        assert len(displacement_y._fc) == 1
        assert displacement_y._fc.time_freq_support.time_frequencies.data == 1
        field = displacement_y._fc[0]
        assert field.component_count == 1
        assert field.data.shape == (121,)
        assert np.allclose(field.data, 0.0)
