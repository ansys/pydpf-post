import ansys.dpf.core as core
from pytest import fixture

import ansys.dpf.post as dpf


@fixture
def simulation(simple_bar):
    return core.Model(simple_bar)


def test_data_frame_from_fields_container(simulation):
    stress_fc = simulation.results.stress().eval()
    df = dpf.DataFrame(data=stress_fc)
    print("\n", df)
