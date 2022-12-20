from ansys.dpf.core.meshed_region import MeshedRegion
import pytest

from ansys.dpf import post


@pytest.fixture()
def simulation(allkindofcomplexity):
    return post.load_solution(allkindofcomplexity, legacy=False)


def test_mesh(simulation):
    assert type(simulation.mesh._meshed_region) == MeshedRegion
    assert len(simulation.mesh._meshed_region.nodes.scoping) >= 15000


def test_available_named_selections(simulation):
    available_named_selections = simulation.mesh.available_named_selections
    assert len(available_named_selections) == 6
    assert "_ELMISC" in available_named_selections
