import pytest

from ansys.dpf import core as dpf
from ansys.dpf.post import connectivity


def test_connectivity_list_iterator():
    property_field = dpf.property_field.PropertyField()
    cli = connectivity.ConnectivityListIterator(conn_list=property_field)
    for i in cli:
        assert isinstance(i, connectivity.ConnectivityListIterator)


def test_connectivity_connectivity_list_idx():
    property_field = dpf.property_field.PropertyField()
    property_field.append(10, scopingid=1)
    property_field.append(11, scopingid=2)
    property_field.append(12, scopingid=3)
    with pytest.raises(ValueError, match="'mode' argument must be"):
        _ = connectivity.ConnectivityListIdx(field=property_field, mode="test")
    # scoping = dpf.mesh_scoping_factory.nodal_scoping([1, 2, 3])
    cli = connectivity.ConnectivityListIdx(
        field=property_field, mode=connectivity.ReturnMode.IDX
    )
    for i in cli:
        assert isinstance(i, list)
    assert cli[1] == [11]
    assert len(cli) == 3
    ref = "ConnectivityListIdx([[10], [11], [12]],__len__=3)"
    assert repr(cli) == ref
    ref = "[[10], [11], [12]]"
    assert str(cli) == ref
    cli2 = cli.by_id
    assert isinstance(cli2, connectivity.ConnectivityListById)

    cli = connectivity.ConnectivityListIdx(
        field=property_field, mode=connectivity.ReturnMode.IDS
    )
    for i in cli:
        assert isinstance(i, list)

    assert [10] in cli


def test_connectivity_connectivity_list_by_id():
    property_field = dpf.property_field.PropertyField()
    property_field.append(10, scopingid=1)
    property_field.append(11, scopingid=2)
    property_field.append(12, scopingid=3)
    cli = connectivity.ConnectivityListById(
        field=property_field,
        mode=connectivity.ReturnMode.IDS,  # , scoping=property_field.scoping
    )
    assert cli[1] == [10]
    for i in cli:
        assert isinstance(i, list)
