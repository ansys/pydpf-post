import ansys.dpf.post as dpf
from ansys.dpf.post.common import AvailableSimulationTypes
from ansys.dpf.post.simulation import (
    HarmonicMechanicalSimulation,
    ModalMechanicalSimulation,
    StaticMechanicalSimulation,
    TransientMechanicalSimulation,
)


def test_load_simulation_static_mechanical(simple_bar, complex_model):
    simulation = dpf.load_simulation(data_sources=simple_bar)
    assert type(simulation) == StaticMechanicalSimulation
    simulation = dpf.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )
    assert type(simulation) == StaticMechanicalSimulation


def test_load_simulation_transient_mechanical(transient_rst, complex_model):
    simulation = dpf.load_simulation(data_sources=transient_rst)
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
