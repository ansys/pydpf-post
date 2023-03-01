import ansys.dpf.core as core
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post.index import (
    CompIndex,
    LabelIndex,
    MeshIndex,
    MultiIndex,
    ResultsIndex,
)
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation


@fixture
def df(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.displacement()


def test_dataframe_core_object(df):
    assert isinstance(df._core_object, core.FieldsContainer)
    assert len(df._core_object) == 1


def test_dataframe_from_fields_container(simple_bar):
    model = core.Model(simple_bar)
    fc = model.results.displacement().eval()
    column_indexes = [
        LabelIndex(name=label, values=fc.get_available_ids_for_label(label))
        for label in fc.labels
    ]
    column_indexes.append(ResultsIndex(values=["U"]))
    column_index = MultiIndex(indexes=column_indexes)

    row_indexes = [
        MeshIndex(location=post.locations.nodal, fc=fc),
        CompIndex(values=["X", "Y", "Z"]),
    ]
    row_index = MultiIndex(indexes=row_indexes)

    df = post.DataFrame(
        data=fc,
        columns=column_index,
        index=row_index,
    )
    assert df.axes == ["node", "comp", "time", "results"]


def test_dataframe_from_error():
    fc = [1, 2, 3]
    with pytest.raises(ValueError, match="not a valid data type"):
        _ = post.DataFrame(data=fc, index=[1, 2])


# def test_dataframe_len(multishells, transient_rst):
#     # Case of one set with two nodal fields
#     model = core.Model(multishells)
#     fc = model.results.stress.on_all_time_freqs.on_location("Nodal").eval()
#     df = post.DataFrame(data=fc)
#     assert len(df) == 2
#     # Case of several sets with one field per set
#     model = core.Model(transient_rst)
#     fc = model.results.displacement.on_all_time_freqs.eval()
#     df = post.DataFrame(data=fc)
#     assert len(df) == 35


def test_dataframe_columns(df):
    columns = df.columns
    print(repr(columns))
    print(columns)


def test_dataframe_index(df):
    index = df.index
    print(repr(index))
    print(index)


def test_dataframe_select_raise(df):
    with pytest.raises(ValueError, match="has no axis"):
        df.select(set_id=1)


def test_dataframe_select(df):
    df2 = df.select(node=[1, 2], time=1, comp="X")
    assert all(df2.mesh_index.values == [1, 2])
    assert df2.index.comp.values == ["X"]
    assert df2.columns.time.values == [1]
    print(df2)


def test_dataframe_iselect(df):
    df2 = df.iselect(node=[0, 1], time=0, comp=0)
    assert all(df2.mesh_index.values == [1, 26])
    assert df2.index.comp.values == ["X"]
    assert df2.columns.time.values == [1]
    print(df2)


def test_dataframe_plot(df):
    df.plot(time=1, node=[1, 2, 3, 4, 5, 6, 7, 8, 9])
