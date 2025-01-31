# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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

    cli = connectivity.ConnectivityListByIndex(
        field=property_field, mode=connectivity.ReturnMode.IDS, scoping=scoping
    )
    ref = "[[100, 101, 102], [102, 103, 104], [104, 105, 106]]"
    assert str(cli) == ref

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
