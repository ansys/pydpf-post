import ansys.dpf.core as core
import pytest
from pytest import fixture

from ansys.dpf import post
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
    df = post.DataFrame(data=fc)
    assert df


def test_dataframe_from_error():
    fc = [1, 2, 3]
    with pytest.raises(ValueError, match="not a valid data type"):
        _ = post.DataFrame(data=fc)


def test_dataframe_len(multishells, transient_rst):
    # Case of one set with two nodal fields
    model = core.Model(multishells)
    fc = model.results.stress.on_all_time_freqs.on_location("Nodal").eval()
    df = post.DataFrame(data=fc)
    assert len(df) == 2
    # Case of several sets with one field per set
    model = core.Model(transient_rst)
    fc = model.results.displacement.on_all_time_freqs.eval()
    df = post.DataFrame(data=fc)
    assert len(df) == 35


def test_dataframe_columns(df):
    columns = df.columns
    print(repr(columns))
    print(columns)


def test_dataframe_index(df):
    index = df.index
    print(repr(index))
    print(index)
