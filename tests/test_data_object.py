import pytest

from ansys.dpf import post


@pytest.fixture()
def displacement(simple_bar):
    return post.load_simulation(simple_bar).displacement(component="Z")


def test_data_object_len(displacement):
    assert len(displacement) == 1


def test_data_object_max(displacement):
    assert displacement.max() == pytest.approx(85061.90584328063)


def test_data_object_max_(displacement):
    assert max(displacement) == pytest.approx(85061.90584328063)


def test_data_object_min(displacement):
    assert displacement.min() == pytest.approx(-1.8264594402081715e-06)


def test_data_object_min_(displacement):
    assert min(displacement) == pytest.approx(-1.8264594402081715e-06)


def test_data_object_as_data_frame(displacement):
    assert False


def test_data_object_as_array(displacement):
    assert False


def test_data_object_plot(displacement, simple_bar):
    displacement.plot(opacity=0.3, shell_layers=None)
    full_displacement = post.load_simulation(simple_bar).displacement()
    displacement.plot(deformation=full_displacement, scale_factor=2.0e4)
    displacement.plot(deformation=displacement, scale_factor=2.0e4)


def test_data_object_animate(displacement):
    assert False
