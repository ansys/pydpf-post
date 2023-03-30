from pytest import fixture

from ansys.dpf import post


@fixture
def fluid_example():
    return r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj"  # noqa


def test_simulation_init(fluid_example):
    fluid_simulation = post.FluidSimulation(fluid_example)
    print(fluid_simulation)
    assert fluid_simulation is not None
