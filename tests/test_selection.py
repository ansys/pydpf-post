import numpy as np

from ansys.dpf import post
from ansys.dpf.post.selection import SpatialSelection


def test_spatial_selection_select_nodes(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection.select_nodes([1, 2, 3])
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.selection.locations.nodal
    assert np.allclose(scoping.ids, [1, 2, 3])

    selection = SpatialSelection()
    selection.select_elements([1, 2, 3, 4])
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.selection.locations.elemental
    assert np.allclose(scoping.ids, [1, 2, 3, 4])
    ids = selection.apply_to(simulation)
    assert np.allclose(ids, [1, 2, 3, 4])


def test_spatial_selection_select_named_selection(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)

    selection = SpatialSelection()
    selection.select_named_selection(
        simulation.mesh.available_named_selections[0],
        location=post.selection.locations.nodal,
    )
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.selection.locations.nodal
    assert scoping.ids.size == 12970
    assert 1857 in scoping.ids
    assert 14826 in scoping.ids
    ids = selection.apply_to(simulation)
    assert len(ids) == 12970
    assert 1857 in ids
    assert 14826 in ids


def test_spatial_selection_intersect(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection1 = SpatialSelection()
    selection1.select_nodes([1, 2, 3])
    _ = selection1._evaluate_on(simulation)

    selection = SpatialSelection()
    selection.select_nodes([1, 2, 3, 4])
    selection.intersect(selection1)
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.selection.locations.nodal
    assert np.allclose(scoping.ids, [1, 2, 3])
    ids = selection.apply_to(simulation)
    assert np.allclose(ids, [1, 2, 3])
