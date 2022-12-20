import ansys.dpf.core as core
from ansys.dpf.post.data_object import DataObject

import numpy as np


def test_data_object_sorting():
    fc = core.FieldsContainer()
    field_1 = core.field_from_array([10, 40, 30, 20])
    field_1.scoping.ids = [1, 4, 3, 2]
    field_2 = core.field_from_array([60, 50, 70, 90, 40])
    field_2.scoping.ids = [6, 5, 7, 9, 4]
    fc.labels = ["mesh_region"]
    fc.add_field({fc.labels[0]: 1}, field_1)
    fc.add_field({fc.labels[0]: 2}, field_2)

    scoping = list(field_1.scoping.ids) + list(field_2.scoping.ids)
    data = DataObject(fc, mesh_scoping=scoping)
    assert data.is_sorted() == False

    data.sort()
    assert data.is_sorted() == True
    assert (data._fc[0].data == np.sort(field_1.data)).all()
    assert (data._fc[1].data == np.sort(field_2.data)).all()
    assert (data._fc[0].scoping.ids == np.sort(field_1.scoping.ids)).all()
    assert (data._fc[1].scoping.ids == np.sort(field_2.scoping.ids)).all()
