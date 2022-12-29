import pytest

from ansys.dpf.post import load_simulation


@pytest.mark.parametrize("example", ["allkindofcomplexity", "simple_bar"])
def test_get_displacements(example, request):
    simulation = load_simulation(request.getfixturevalue(example))
    # simulation = load_simulation(example)

    for field in simulation.results:
        if field.name == "displacement":
            n_dim = field.n_components
            break

    n_time_steps = len(simulation.time_freq_support.time_frequencies.data)
    # n_nodes = simulation.mesh._meshed_region.nodes.n_nodes

    disp = simulation.displacement()
    assert len(disp) == n_time_steps
    assert disp[0].n_dim == n_dim
    # assert disp[0].elementary_data_count == n_nodes # HOW CAN WE CHECK THE n_entities?
    disp_x = simulation.displacement(component=0)
    assert len(disp_x) == n_time_steps
    assert disp_x[0].n_dim == 1
    # assert disp[0].elementary_data_count == n_nodes # HOW CAN WE CHECK THE n_entities?

    if n_time_steps > 1:
        disp_time = simulation.displacement(steps=[1, 25])
        assert len(disp_time) == 2
        assert disp[0].n_dim == n_dim
        # assert disp_time[0].elementary_data_count == n_nodes # HOW CAN WE CHECK THE n_entities?


def test_get_stresses(allkindofcomplexity, static_rst):
    ###########################################################################
    # allkindofcomplexity example
    simulation = load_simulation(allkindofcomplexity)

    for field in simulation.results:
        if field.name == "stress":
            n_dim = field.n_components
            break

    # n_elems = simulation.mesh._meshed_region.elements.n_elements
    # n_nodes = simulation.mesh._meshed_region.nodes.n_nodes
    n_stress_bodies = 2
    # Nodal stress
    stress = simulation.nodal_stress()
    assert len(stress) == n_stress_bodies
    assert stress[0].n_dim == n_dim

    stress_xy_nodes = simulation.nodal_stress(nodes=[100, 101, 102, 103])
    assert len(stress_xy_nodes) == 1
    assert stress_xy_nodes[0].n_dim == n_dim
    assert stress_xy_nodes[0].n_data == 4

    # Elemental stress
    stress = simulation.elemental_stress()
    assert len(stress) == n_stress_bodies
    assert stress[0].n_dim == n_dim

    # Raw stress
    stress = simulation.raw_stress()
    assert len(stress) == 1
    # WHERE CAN WE GET BETTER INFO ON WHAT'S ON EACH ENTRY OF THE FIELDSCONTAINER?
    assert stress[0].n_dim == n_dim

    ###########################################################################
    # static_rst example
    simulation = load_simulation(static_rst)

    for field in simulation.results:
        if field.name == "stress":
            n_dim = field.n_components
            break

    stress = simulation.nodal_stress()
    assert len(stress) == 1
    assert stress[0].n_dim == n_dim

    stress = simulation.elemental_stress()
    assert len(stress) == 1
    assert stress[0].n_dim == n_dim


components = (
    "X",
    "XY",
    "XZ",
    "Y",
    "YZ",
    "Z",
    "test",
    pytest.param(0.1, marks=pytest.mark.xfail(strict=True, raises=TypeError)),
)


@pytest.mark.parametrize("components", components)
def test_get_stresses_by_component(components, allkindofcomplexity):
    simulation = load_simulation(allkindofcomplexity)
    for field in simulation.results:
        if field.name == "stress":
            n_dim = field.n_components
            break

    for component in components:
        stress = simulation.nodal_stress(component=component)
        assert stress[0].n_dim == n_dim if component == "test" else 1


def test_data_frame():
    from ansys.dpf.post import examples, load_simulation

    simulation = load_simulation(examples.download_transient_result())
    dispObject = simulation.displacement(nodes=[20, 200, 400], steps=[25])
    df = dispObject.as_data_frame()
    assert len(df.columns) == 3
    assert len(df.index) == 3

    stressObject = simulation.nodal_stress()
    df = stressObject.as_data_frame()
    assert len(df.columns) == 35 * 6
    assert len(df.index) == 3812
