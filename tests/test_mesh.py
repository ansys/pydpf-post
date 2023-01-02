from ansys.dpf.core.meshed_region import MeshedRegion
import pytest

from ansys.dpf import post


@pytest.fixture()
def simulation(allkindofcomplexity):
    return post.load_simulation(allkindofcomplexity)


def test_mesh(simulation):
    mesh = simulation.mesh
    assert type(mesh._meshed_region) == MeshedRegion
    assert len(mesh._meshed_region.nodes.scoping) >= 15000

    available_named_selections = simulation.mesh.available_named_selections
    assert len(available_named_selections) == 6
    assert "_ELMISC" in available_named_selections

    available_property_fields = mesh.available_property_fields
    assert len(available_property_fields) == 6
    assert "connectivity" in available_property_fields

    assert mesh.grid.n_cells == 10292
    assert mesh.grid.n_points == 15129

    assert mesh.elements.n_elements == 10292
    assert mesh.nodes.n_nodes == 15129
    assert mesh.unit == "m"


def test_mesh_plot(simulation):
    mesh = simulation.mesh
    displacement = simulation.displacement()
    # mesh.plot()
    # displacement.plot()
    print(displacement)
    mesh.plot(data=displacement)
    mesh.plot(deformation=displacement, scale_factor=2.0)
