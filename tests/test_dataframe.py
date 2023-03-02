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
from ansys.dpf.post.transient_mechanical_simulation import TransientMechanicalSimulation


@fixture
def df(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.displacement()


@fixture
def elastic_strain_df(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.elastic_strain_nodal()


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
    # print(df)
    df2 = df.select(node=[1, 2], time=1, comp="X")
    assert all(df2.mesh_index.values == [1, 2])
    assert df2.index.comp.values == ["X"]
    assert df2.columns.time.values == [1]
    # print(df2)


def test_dataframe_iselect(df):
    df2 = df.iselect(node=[0, 1], time=0, comp=0)
    assert all(df2.mesh_index.values == [1, 26])
    assert df2.index.comp.values == ["X"]
    assert df2.columns.time.values == [1]
    # print(df2)


def test_dataframe_plot(df):
    df.plot(time=1, node=[1, 2, 3, 4, 5, 6, 7, 8, 9])


def test_dataframe_repr(df):
    ref = (
        "DataFrame<index=MultiIndex<[MeshIndex<name=\"node\", dtype=<class 'int'>>, "
        "Index<name=\"comp\", dtype=<class 'str'>>]>, columns=MultiIndex<[ResultIndex<['U']>, "
        "Index<name=\"time\", dtype=<class 'int'>>]>>"
    )
    assert repr(df) == ref


def test_dataframe_str(transient_rst):
    simulation = TransientMechanicalSimulation(transient_rst)
    df = simulation.displacement(all_sets=True)
    # print(df)
    ref = """  
             results         U                                                  
                time         1         2         3         4         5         6
      node      comp                                                            
       525         X  0.00e+00  4.85e-05  2.30e-04  6.51e-04  1.48e-03  2.93e-03
                   Y  0.00e+00  2.87e-04  1.14e-03  2.54e-03  4.41e-03  6.59e-03
                   Z  0.00e+00 -1.26e-10 -4.34e-10 -8.29e-10 -1.15e-09 -1.39e-09
       534         X  0.00e+00  6.55e-06  1.05e-04  5.30e-04  1.67e-03  4.02e-03
                   Y  0.00e+00  6.27e-04  2.51e-03  5.62e-03  9.86e-03  1.50e-02
"""  # noqa: W291
    assert str(df) == ref


def test_dataframe_str_comp(df):
    # 3D str
    stri = str(df)
    expected_strs = ["X", "Y", "Z", "results", "U", "node"]
    for expected_str in expected_strs:
        assert expected_str in stri
    # 1D str
    stri = str(df.select(comp="X"))
    expected_strs = ["X", "results", "U", "node"]
    for expected_str in expected_strs:
        assert expected_str in stri
    assert "Y" not in stri


def test_dataframe_str_tensor(elastic_strain_df):
    # 3D str
    stri = str(elastic_strain_df)
    expected_strs = ["XX", "XY", "XZ", "results", "EPEL", "node"]
    for expected_str in expected_strs:
        assert expected_str in stri
    # 1D str
    stri = str(elastic_strain_df.select(comp=["XX", "XY"]))
    expected_strs = ["XX", "XY", "results", "U", "node"]
    for expected_str in expected_strs:
        assert expected_str in stri
    assert "ZZ" not in stri