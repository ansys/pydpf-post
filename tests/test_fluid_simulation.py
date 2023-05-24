from ansys.dpf.core import examples
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post


@fixture
def fluid_simulation():
    fluid_example_files = examples.download_fluent_files()
    ds = dpf.DataSources()
    ds.set_result_file_path(
        fluid_example_files["cas"],
        key="cas",
    )
    ds.add_file_path(
        fluid_example_files["dat"],
        key="dat",
    )
    return post.FluidSimulation(ds)  # noqa


def test_simulation_init(fluid_simulation):
    assert fluid_simulation is not None


def test_density(fluid_simulation):
    result = fluid_simulation.density()
    assert result is not None
    print(result)
