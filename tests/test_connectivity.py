import pytest

from ansys.dpf import core as dpf
from ansys.dpf.post import connectivity


def test_connectivity_connectivity_list_idx():
    property_field = dpf.property_field.PropertyField()
    property_field.append([0, 1, 2], scopingid=1)
    property_field.append([2, 3, 4], scopingid=2)
    property_field.append([4, 5, 6], scopingid=3)
    scoping = dpf.mesh_scoping_factory.nodal_scoping(
        [100, 101, 102, 103, 104, 105, 106]
    )
    with pytest.raises(ValueError, match="'mode' argument must be"):
        _ = connectivity.ConnectivityListByIndex(
            field=property_field, mode="test", scoping=scoping
        )
    cli = connectivity.ConnectivityListByIndex(
        field=property_field, mode=connectivity.ReturnMode.IDX, scoping=scoping
    )
    for i in cli:
        assert isinstance(i, list)
    assert cli[1] == [2, 3, 4]
    assert len(cli) == 3
    ref = "ConnectivityListByIndex([[0 1 2], [2 3 4], [4 5 6]], __len__=3)"
    assert repr(cli) == ref
    ref = "[[0 1 2], [2 3 4], [4 5 6]]"
    assert str(cli) == ref
    cli2 = cli.by_id
    assert isinstance(cli2, connectivity.ConnectivityListById)

    property_field.append([6, 7, 8], scopingid=4)
    scoping = dpf.mesh_scoping_factory.nodal_scoping(
        [100, 101, 102, 103, 104, 105, 106, 107, 108]
    )
    cli = connectivity.ConnectivityListByIndex(
        field=property_field, mode=connectivity.ReturnMode.IDS, scoping=scoping
    )
    for i in cli:
        assert isinstance(i, list)
    for i in cli:
        assert isinstance(i, list)
    ref = "ConnectivityListByIndex([[100, 101, 102], ..., [106, 107, 108]], __len__=4)"
    assert repr(cli) == ref
    ref = "[[100, 101, 102], ..., [106, 107, 108]]"
    assert str(cli) == ref


def test_connectivity_connectivity_list_by_id():
    property_field = dpf.property_field.PropertyField()
    property_field.append([0, 1, 2], scopingid=1)
    property_field.append([2, 3, 4], scopingid=2)
    property_field.append([4, 5, 6], scopingid=3)
    property_field.append([6, 7, 8], scopingid=4)
    scoping = dpf.mesh_scoping_factory.nodal_scoping(
        [100, 101, 102, 103, 104, 105, 106, 107, 108]
    )
    cli = connectivity.ConnectivityListById(
        field=property_field, mode=connectivity.ReturnMode.IDS, scoping=scoping
    )
    assert cli[1] == [100, 101, 102]
    for i in cli:
        assert isinstance(i, list)
    for i in cli:
        assert isinstance(i, list)
    assert len(cli) == 4
    ref = "ConnectivityListById([[100, 101, 102], ..., [106, 107, 108]], __len__=4)"
    assert repr(cli) == ref
    ref = "[[100, 101, 102], ..., [106, 107, 108]]"
    assert str(cli) == ref

    cli = connectivity.ConnectivityListById(
        field=property_field, mode=connectivity.ReturnMode.IDS, scoping=scoping
    )
    for i in cli:
        assert isinstance(i, list)
    for i in cli:
        assert isinstance(i, list)
