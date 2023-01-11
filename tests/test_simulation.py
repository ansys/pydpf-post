from pytest import fixture

import ansys.dpf.post as dpf


@fixture
def simulation(simple_bar):
    return dpf.load_simulation(data_sources=simple_bar)


def test_simulation_results(simulation):
    results = simulation.results
    assert len(results) == 10
    assert all(isinstance(x, str) for x in results)


def test_simulation_geometries(simulation):
    geometries = simulation.geometries
    assert geometries == []


def test_simulation_boundary_conditions(simulation):
    boundary_conditions = simulation.boundary_conditions
    assert boundary_conditions == []


def test_simulation_loads(simulation):
    loads = simulation.loads
    assert loads == []


def test_simulation_mesh(simulation):
    mesh = simulation.mesh
    assert isinstance(mesh, dpf.mesh.Mesh)


def test_simulation_named_selections(simulation):
    named_selections = simulation.named_selections
    assert len(named_selections) == 1
    assert all(isinstance(x, str) for x in named_selections)


def test_simulation_active_selection(simulation):
    assert simulation.active_selection is None
    selection = dpf.selection.Selection()
    simulation.activate_selection(selection=selection)
    assert simulation.active_selection == selection
    simulation.deactivate_selection()
    assert simulation.active_selection is None
