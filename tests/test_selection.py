import numpy as np
from ansys.dpf.core.common import locations

from ansys.dpf import post
from ansys.dpf.post.selection import Selection, SpatialSelection
import numpy


def test_spatial_scoping_selection(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    selection = SpatialSelection()
    selection.select_nodes([1, 2, 3])
    scoping = selection.evaluate_on(solution)
    assert scoping.location == post.selection.locations.nodal
    assert np.allclose(scoping.ids, [1, 2, 3])

    selection = SpatialSelection()
    selection.select_elements([1, 2, 3, 4])
    scoping = selection.evaluate_on(solution)
    assert scoping.location == post.selection.locations.elemental
    assert np.allclose(scoping.ids, [1, 2, 3, 4])


def test_spatial_intersect_selection(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    selection1 = SpatialSelection()
    selection1.select_nodes([1, 2, 3])
    scoping = selection1.evaluate_on(solution)

    selection = SpatialSelection()
    selection.select_nodes([1, 2, 3, 4])
    selection.select_intersection_with(selection1)
    scoping = selection.evaluate_on(solution)
    assert scoping.location == post.selection.locations.nodal
    assert np.allclose(scoping.ids, [1, 2, 3])
