import ansys.dpf.core as core
from pytest import fixture

from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation


@fixture
def df(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.displacement()


def test_dataframe_core_object(df):
    assert isinstance(df._core_object, core.FieldsContainer)
    assert len(df._core_object) == 1
