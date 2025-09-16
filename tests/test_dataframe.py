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

import ansys.dpf.core as core
import numpy as np
import pytest
from pytest import fixture

from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.index import (
    CompIndex,
    LabelIndex,
    MeshIndex,
    MultiIndex,
    ResultsIndex,
    ref_labels,
)
from ansys.dpf.post.modal_mechanical_simulation import ModalMechanicalSimulation
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation
from ansys.dpf.post.transient_mechanical_simulation import TransientMechanicalSimulation
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0


@fixture
def df(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    # print(simulation._model._server.version)
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
    repr(columns)
    # print(repr(columns))
    str(columns)
    # print(columns)


def test_dataframe_index(df):
    index = df.index
    # print(repr(index))
    repr(index)
    # print(index)
    str(index)


def test_dataframe_select_raise(df, transient_rst):
    # with pytest.raises(ValueError, match="has no axis"):
    #     df.select(toto=1)

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


def test_dataframe_select_no_set_index():
    simulation = post.StaticMechanicalSimulation(examples.find_simple_bar())
    df = simulation.mesh.coordinates
    df2 = df.select(components="X")
    assert len(df2.index) == 2
    assert df2.index.components.values == ["X"]
    assert len(df2.mesh_index.values) == len(df.mesh_index.values)


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
def test_dataframe_select_cells():
    simulation = post.FluidSimulation(examples.fluid_axial_model())
    df = simulation.enthalpy()
    df.select(cell_ids=[1])


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
def test_dataframe_select_with_labels():
    fluid_file = examples.download_cfx_mixing_elbow()
    simulation = post.FluidSimulation(fluid_file)
    df = simulation.enthalpy()
    df2 = df.select(node_ids=[1])
    ref = """
  results        H_S (J/kg)
  set_ids                 1
    phase Water at 25 C (2)
 node_ids                  
        1       -1.5921e+07
"""  # noqa: W291, E501
    assert str(df2) == ref


def test_dataframe_iselect(df):
    df2 = df.iselect(node_ids=[0, 1], set_ids=[0], components=0)
    assert all(df2.mesh_index.values == [1, 26])
    assert df2.index.components.values == ["X"]
    assert df2.columns.set_ids.values == [1]
    # print(df2)


# def test_dataframe_plot(df, multi_stage_cyclic):
#     df.plot(set_ids=1, node_ids=[1, 2, 3, 4, 5, 6, 7, 8, 9])

#     simulation = post.ModalMechanicalSimulation(multi_stage_cyclic)
#     # df2 = simulation.displacement(expand_cyclic=False)  # TODO fix plot bug
#     df2 = simulation.stress_nodal(expand_cyclic=False)
#     # print(df2)
#     df2.plot()
#     df2.plot(stage=0)
#     with pytest.raises(ValueError, match="must be a single value"):
#         df2.plot(stage=[0, 1])


def test_dataframe_plot_warn(df):
    with pytest.warns(UserWarning, match="did not return data"):
        plt = df.plot(set_ids=99)
        assert plt is None


# def test_dataframe_animate(transient_rst):
#     simulation = TransientMechanicalSimulation(transient_rst)
#     # Animate displacement
#     df = simulation.displacement(all_sets=True)
#     df.animate()
#     df.animate(scale_factor=5.0, deform=True)

#     df = simulation.stress_nodal(all_sets=True)
#     df.animate(deform=True, scale_factor=5.0)


def test_dataframe_repr(df):
    ref = (
        "DataFrame<"
        "index=MultiIndex<["
        "MeshIndex<name=\"node_ids\", dtype=<class 'int'>>, "
        "CompIndex<name=\"components\", dtype=<class 'str'>>"
        "]>, "
        "columns=MultiIndex<["
        "ResultIndex<['U (m)']>, "
        "SetIndex<values=[1]>"
        "]>>"
    )
    assert repr(df) == ref


def test_dataframe_str(transient_rst, modal_frame):
    simulation = TransientMechanicalSimulation(transient_rst)
    df = simulation.displacement(all_sets=True)
    # print(df)
    ref = """
             results      U (m)                                                             ...
             set_ids          1           2           3           4           5           6 ...
 node_ids components                                                                        ...
      525          X 0.0000e+00  4.8506e-05  2.3022e-04  6.5140e-04  1.4812e-03  2.9324e-03 ...
                   Y 0.0000e+00  2.8732e-04  1.1437e-03  2.5408e-03  4.4069e-03  6.5936e-03 ...
                   Z 0.0000e+00 -1.2615e-10 -4.3450e-10 -8.2924e-10 -1.1459e-09 -1.3910e-09 ...
      534          X 0.0000e+00  6.5467e-06  1.0495e-04  5.3050e-04  1.6666e-03  4.0153e-03 ...
                   Y 0.0000e+00  6.2670e-04  2.5072e-03  5.6168e-03  9.8601e-03  1.4993e-02 ...
                   Z 0.0000e+00 -3.1963e-10 -1.1039e-09 -2.1288e-09 -2.9359e-09 -3.5675e-09 ...
      ...        ...        ...         ...         ...         ...         ...         ... ...
"""  # noqa: W291, E501
    assert str(df) == ref
    df2 = df.select(node_ids=525)
    # print(df2)
    ref = """
             results      U (m)                                                             ...
             set_ids          1           2           3           4           5           6 ...
 node_ids components                                                                        ...
      525          X 0.0000e+00  4.8506e-05  2.3022e-04  6.5140e-04  1.4812e-03  2.9324e-03 ...
                   Y 0.0000e+00  2.8732e-04  1.1437e-03  2.5408e-03  4.4069e-03  6.5936e-03 ...
                   Z 0.0000e+00 -1.2615e-10 -4.3450e-10 -8.2924e-10 -1.1459e-09 -1.3910e-09 ...
"""  # noqa: W291, E501
    assert str(df2) == ref

    df = simulation.stress()
    # print(df)
    # print(df._fc[0].get_entity_data_by_id(391))
    ref = """
                results      S (Pa)                                                             ...
                set_ids          35                                                             ...
                   node           0           1           2           3           4           5 ...
 element_ids components                                                                         ...
         391         XX -3.2780e+05 -4.6382e+06 -2.3568e+07 -3.9276e+07  1.6355e+07  1.8076e+07 ...
                     YY  1.3601e+06  1.2931e+06 -3.5347e+07 -2.7237e+07  2.8319e+07  5.2558e+06 ...
                     ZZ  1.4909e+08  1.2041e+08  2.0150e+08  1.8145e+08  1.0508e+08  7.7621e+07 ...
                     XY -4.8869e+06 -6.0662e+06 -5.2336e+06 -3.7544e+06 -1.3022e+07 -7.5306e+06 ...
                     YZ  1.4304e+07  2.3483e+07 -2.8879e+07 -3.6248e+06  5.1991e+05  2.4472e+06 ...
                     XZ  1.6546e+07  1.7723e+07 -6.1648e+06 -3.2608e+07  8.8243e+06  2.9268e+06 ...
         ...        ...         ...         ...         ...         ...         ...         ... ...
"""  # noqa: W291, E501
    assert str(df) == ref
    # print(df._fc[0].get_entity_data_by_id(391)[0][1])
    # print(df._fc[0].get_entity_data_by_id(456)[0][1])
    # print(df._fc[0].get_entity_data_by_id(326)[0][1])
    # print(df._fc[0].get_entity_data_by_id(586)[0][1])
    # print(df._fc[0].get_entity_data_by_id(261)[0][1])
    # print(df._fc[0].get_entity_data_by_id(66)[0][1])
    df2 = df.select(components="YY")
    # print(df2)
    ref = """
                results      S (Pa)                                                             ...
                set_ids          35                                                             ...
                   node           0           1           2           3           4           5 ...
 element_ids components                                                                         ...
         391         YY  1.3601e+06  1.2931e+06 -3.5347e+07 -2.7237e+07  2.8319e+07  5.2558e+06 ...
         456             1.9918e+05 -8.0343e+06 -1.6874e+07 -3.3505e+07  5.8703e+06  1.9095e+06 ...
         326             1.9563e+06  1.8186e+06 -2.6996e+07 -1.2659e+07  1.8008e+07  2.7146e+07 ...
         586            -5.3043e+06  7.5858e+06 -1.7397e+07  5.2490e+05 -1.0581e+06 -4.3758e+07 ...
         261            -2.4112e+05 -1.3009e+06  2.4026e+05  6.0941e+05 -5.3125e+06  4.1034e+07 ...
          66             4.8877e+05 -9.3109e+05 -1.5726e+06 -3.0581e+05 -6.1510e+06 -3.1282e+07 ...
         ...        ...         ...         ...         ...         ...         ...         ... ...
"""  # noqa: W291, E501
    assert str(df2) == ref

    df = simulation.stress(all_sets=True)
    # print(df)
    ref = """
                results     S (Pa)                                                        ...
                set_ids          1                                                        ...
                   node          0          1          2          3          4          5 ...
 element_ids components                                                                   ...
         391         XX 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 ...
                     YY 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 ...
                     ZZ 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 ...
                     XY 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 ...
                     YZ 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 ...
                     XZ 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 ...
         ...        ...        ...        ...        ...        ...        ...        ... ...
"""  # noqa: W291, E501
    assert str(df) == ref
    assert df.display_max_columns == 6
    assert df.display_max_rows == 6
    df.display_max_columns = 2
    df.display_max_rows = 2
    ref = """
                results     S (Pa)            ...
                set_ids          1            ...
                   node          0          1 ...
 element_ids components                       ...
         391         XX 0.0000e+00 0.0000e+00 ...
                     YY 0.0000e+00 0.0000e+00 ...
         ...        ...        ...        ... ...
"""  # noqa: W291, E501
    assert str(df) == ref
    # Test for cyclic with base sector label
    modal_simulation = ModalMechanicalSimulation(examples.find_simple_cyclic())
    df = modal_simulation.displacement(expand_cyclic=False)
    ref = """
              results       U (m)
              set_ids           1
          base_sector           1
 node_ids  components            
        1           X  4.9812e-13
                    Y  2.4100e+02
                    Z  8.9709e-12
       14           X -1.9511e-12
                    Y  1.9261e+02
                    Z  5.0359e-12
      ...         ...         ...
"""  # noqa: W291, E501
    assert str(df) == ref


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


def test_dataframe_array(elastic_strain_df):
    arr = elastic_strain_df.array
    assert isinstance(arr, np.ndarray)
    assert arr.ndim == 2
    assert arr.size == 486
    assert len(arr[0]) == 6
    assert arr[0][0] - 1.5236e-07 < 1.0e-13


def test_dataframe_array_raise(transient_rst):
    simulation = TransientMechanicalSimulation(transient_rst)
    df = simulation.displacement(all_sets=True)
    with pytest.raises(
        ValueError, match="Can only export to array if the DataFrame contains a single"
    ):
        _ = df.array


def test_dataframe_min_max():
    simulation = post.TransientMechanicalSimulation(examples.download_crankshaft())
    df = simulation.displacement(all_sets=True)
    # Over the mesh entities
    min_over_mesh = [[-0.00074732, -0.00040138, -0.00021555]]
    assert np.all(np.isclose(df.min()._fc[0].data.tolist(), min_over_mesh))
    assert np.all(np.isclose(df.min(axis=0)._fc[0].data.tolist(), min_over_mesh))
    assert np.all(
        np.isclose(df.min(axis="node_ids")._fc[0].data.tolist(), min_over_mesh)
    )

    max_over_mesh = [[0.00073303, 0.00139618, 0.00021567]]
    assert np.all(np.isclose(df.max()._fc[0].data.tolist(), max_over_mesh))
    assert np.all(np.isclose(df.max(axis=0)._fc[0].data.tolist(), max_over_mesh))
    assert np.all(
        np.isclose(df.max(axis="node_ids")._fc[0].data.tolist(), max_over_mesh)
    )

    # Over the SetIndex
    min_over_time = [-3.41368775e-05, 5.16665595e-04, -4.13456506e-06]
    assert np.all(np.isclose(df.min(axis=1)._fc[0].data[0].tolist(), min_over_time))
    assert np.all(
        np.isclose(df.min(axis="set_ids")._fc[0].data[0].tolist(), min_over_time)
    )
    max_over_time = [5.67807472e-06, 1.54174694e-03, -2.63976203e-06]
    assert np.all(np.isclose(df.max(axis=1)._fc[0].data[0].tolist(), max_over_time))
    assert np.all(
        np.isclose(df.max(axis="set_ids")._fc[0].data[0].tolist(), max_over_time)
    )

    # Raise unrecognized axis
    with pytest.raises(ValueError, match="is not an available axis value"):
        df.max(axis="raises")
