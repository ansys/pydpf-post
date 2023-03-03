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
    ref_labels,
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
    assert df.index.names == [ref_labels.node_ids, ref_labels.components]
    assert df.columns.names == [ref_labels.time, ref_labels.results]


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


def test_dataframe_select_raise(df, transient_rst):
    with pytest.raises(ValueError, match="has no axis"):
        df.select(toto=1)

    with pytest.raises(NotImplementedError, match="Element selection"):
        simulation = TransientMechanicalSimulation(transient_rst)
        df = simulation.stress()
        _ = df.select(element_ids=391)


def test_dataframe_select(df):
    # print(df)
    df2 = df.select(node_ids=[1, 2], set_ids=1, components="X")
    assert all(df2.mesh_index.values == [1, 2])
    assert df2.index.components.values == ["X"]
    assert df2.columns.set_ids.values == [1]
    # print(df2)


def test_dataframe_iselect(df):
    df2 = df.iselect(node_ids=[0, 1], set_ids=[0], components=0)
    assert all(df2.mesh_index.values == [1, 26])
    assert df2.index.components.values == ["X"]
    assert df2.columns.set_ids.values == [1]
    print(df2)


def test_dataframe_plot(df):
    df.plot(set_ids=1, node_ids=[1, 2, 3, 4, 5, 6, 7, 8, 9])


def test_dataframe_animate(transient_rst):
    simulation = TransientMechanicalSimulation(transient_rst)
    # Animate displacement
    df = simulation.displacement(all_sets=True)
    # df.animate()
    df.animate(scale_factor=5.0, deform=True, save_as="test_dataframe_animate.gif")
    # Animate nodal stress -> Does not work
    df2 = simulation.stress_nodal(all_sets=True)
    # df2.animate()
    # assert False


def test_dataframe_repr(df):
    ref = (
        "DataFrame<index=MultiIndex<[MeshIndex<name=\"node_ids\", dtype=<class 'int'>>, "
        "Index<name=\"components\", dtype=<class 'str'>>]>, columns=MultiIndex<[ResultIndex<['U']>, "  # noqa: E501
        "SetIndex<values=[1]>]>>"
    )
    assert repr(df) == ref


def test_dataframe_str(transient_rst):
    simulation = TransientMechanicalSimulation(transient_rst)
    df = simulation.displacement(all_sets=True)
    print(df)
    ref = """
                 results           U                                                                     ...
                 set_ids           1           2           3           4           5           6         ...
    node_ids  components                                                                                 ...
         525           X  0.0000e+00  4.8506e-05  2.3022e-04  6.5140e-04  1.4812e-03  2.9324e-03         ...
                       Y  0.0000e+00  2.8732e-04  1.1437e-03  2.5408e-03  4.4069e-03  6.5936e-03         ...
                       Z  0.0000e+00 -1.2615e-10 -4.3450e-10 -8.2924e-10 -1.1459e-09 -1.3910e-09         ...
         534           X  0.0000e+00  6.5467e-06  1.0495e-04  5.3050e-04  1.6666e-03  4.0153e-03         ...
                       Y  0.0000e+00  6.2670e-04  2.5072e-03  5.6168e-03  9.8601e-03  1.4993e-02         ...
                       Z  0.0000e+00 -3.1963e-10 -1.1039e-09 -2.1288e-09 -2.9359e-09 -3.5675e-09         ...
         ...
"""  # noqa: W291, E501
    assert str(df) == ref
    df2 = df.select(node_ids=525)
    print(df2)
    ref = """
                 results           U                                                                     ...
                 set_ids           1           2           3           4           5           6         ...
    node_ids  components                                                                                 ...
         525           X  0.0000e+00  4.8506e-05  2.3022e-04  6.5140e-04  1.4812e-03  2.9324e-03         ...
                       Y  0.0000e+00  2.8732e-04  1.1437e-03  2.5408e-03  4.4069e-03  6.5936e-03         ...
                       Z  0.0000e+00 -1.2615e-10 -4.3450e-10 -8.2924e-10 -1.1459e-09 -1.3910e-09         ...
"""  # noqa: W291, E501
    assert str(df2) == ref

    df = simulation.stress()
    print(df)
    print(df._fc[0].get_entity_data_by_id(391))
    ref = """
                 results           S
                 set_ids          35
 element_ids  components            
         391      XX (1) -3.2780e+05
                  YY (1)  1.3601e+06
                  ZZ (1)  1.4909e+08
                  XY (1) -4.8869e+06
                  YZ (1)  1.4304e+07
                  XZ (1)  1.6546e+07
         ...
"""  # noqa: W291, E501
    assert str(df) == ref
    df2 = df.select(components="YY")
    print(df2)
    ref = """
                 results           S
                 set_ids          35
 element_ids  components            
         391      YY (1)  1.3601e+06
         391      YY (2)  1.2931e+06
         391      YY (3) -3.5347e+07
         391      YY (4) -2.7237e+07
         391      YY (5)  2.8319e+07
         391      YY (6)  5.2558e+06
         ...
"""  # noqa: W291, E501
    assert str(df2) == ref


def test_dataframe_str_comp(df):
    # 3D str
    stri = str(df)
    expected_strs = ["X", "Y", "Z", ref_labels.results, "U", ref_labels.node_ids]
    for expected_str in expected_strs:
        assert expected_str in stri
    # 1D str
    stri = str(df.select(components="X"))
    expected_strs = ["X", ref_labels.results, "U", ref_labels.node_ids]
    for expected_str in expected_strs:
        assert expected_str in stri
    assert "Y" not in stri


def test_dataframe_str_tensor(elastic_strain_df):
    # 3D str
    stri = str(elastic_strain_df)
    expected_strs = ["XX", "XY", "XZ", ref_labels.results, "EPEL", ref_labels.node_ids]
    for expected_str in expected_strs:
        assert expected_str in stri
    # 1D str
    stri = str(elastic_strain_df.select(components=["XX", "XY"]))
    expected_strs = ["XX", "XY", ref_labels.results, "EPEL", ref_labels.node_ids]
    for expected_str in expected_strs:
        assert expected_str in stri
    assert "ZZ" not in stri
