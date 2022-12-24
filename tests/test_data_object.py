import numpy as np

from ansys.dpf import post


def test_get_data_by_id(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    disp = simulation.displacement()
    # Request all IDs
    data, ids = disp.get_data_by_id()
    assert len(data) == len(ids)
    # Request a few ID numbers
    data, ids = disp.get_data_by_id(get_ids=[1, 50, 100])
    assert len(data) == len(ids)
    assert (ids == np.array([1, 50, 100])).all()
    # Request some IDs that are not available
    data, ids = disp.get_data_by_id(get_ids=[1, 50, 100, -1])
    assert len(data) == len(ids)
    assert (ids == np.array([1, 50, 100])).all()
