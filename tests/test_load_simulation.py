import pytest

import ansys.dpf.post as dpf
from ansys.dpf.post.common import AvailableSimulationTypes
from ansys.dpf.post.modal_mechanical_simulation import ModalMechanicalSimulation
from ansys.dpf.post.simulation import HarmonicMechanicalSimulation
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation
from ansys.dpf.post.transient_mechanical_simulation import TransientMechanicalSimulation


def test_load_simulation_static_mechanical(simple_bar, complex_model):
    simulation = dpf.load_simulation(data_sources=simple_bar)
    assert type(simulation) == StaticMechanicalSimulation
    simulation = dpf.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )
    assert type(simulation) == StaticMechanicalSimulation


def test_load_simulation_transient_mechanical(plate_msup, complex_model):
    simulation = dpf.load_simulation(data_sources=plate_msup)
    assert type(simulation) == TransientMechanicalSimulation
    simulation = dpf.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )
    assert type(simulation) == TransientMechanicalSimulation


def test_load_simulation_modal_mechanical(modalallkindofcomplexity, complex_model):
    simulation = dpf.load_simulation(data_sources=modalallkindofcomplexity)
    assert type(simulation) == ModalMechanicalSimulation
    simulation = dpf.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )
    assert type(simulation) == TransientMechanicalSimulation


def test_load_simulation_harmonic_mechanical(complex_model, simple_bar):
    simulation = dpf.load_simulation(data_sources=complex_model)
    assert type(simulation) == HarmonicMechanicalSimulation
    simulation = dpf.load_simulation(
        data_sources=simple_bar,
        simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    )
    assert type(simulation) == HarmonicMechanicalSimulation


def test_load_simulation_raise_simulation_type(simple_bar):
    with pytest.raises(ValueError, match="is not a recognized simulation type"):
        _ = dpf.load_simulation(
            data_sources=simple_bar,
            simulation_type="test",
        )
