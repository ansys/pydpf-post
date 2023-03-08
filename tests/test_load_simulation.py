import pytest

from ansys.dpf import post
from ansys.dpf.post.common import AvailableSimulationTypes
from ansys.dpf.post.harmonic_mechanical_simulation import HarmonicMechanicalSimulation
from ansys.dpf.post.modal_mechanical_simulation import ModalMechanicalSimulation
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation
from ansys.dpf.post.transient_mechanical_simulation import TransientMechanicalSimulation


def test_load_simulation_static_mechanical(simple_bar, complex_model):
    simulation = post.load_simulation(data_sources=simple_bar)
    assert type(simulation) == StaticMechanicalSimulation
    simulation = post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )
    assert isinstance(simulation, StaticMechanicalSimulation)


def test_load_simulation_transient_mechanical(plate_msup, complex_model):
    simulation = post.load_simulation(data_sources=plate_msup)
    assert type(simulation) == TransientMechanicalSimulation
    simulation = post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )
    assert isinstance(simulation, TransientMechanicalSimulation)


def test_load_simulation_modal_mechanical(modalallkindofcomplexity, complex_model):
    simulation = post.load_simulation(data_sources=modalallkindofcomplexity)
    assert type(simulation) == ModalMechanicalSimulation
    simulation = post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )
    assert isinstance(simulation, TransientMechanicalSimulation)


def test_load_simulation_harmonic_mechanical(complex_model, simple_bar):
    simulation = post.load_simulation(data_sources=complex_model)
    assert type(simulation) == HarmonicMechanicalSimulation
    simulation = post.load_simulation(
        data_sources=simple_bar,
        simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    )
    assert isinstance(simulation, HarmonicMechanicalSimulation)


def test_load_simulation_raise_simulation_type(simple_bar):
    with pytest.raises(ValueError, match="is not a recognized simulation type"):
        _ = post.load_simulation(
            data_sources=simple_bar,
            simulation_type="test",
        )
