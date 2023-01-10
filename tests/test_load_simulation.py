import ansys.dpf.post as dpf
from ansys.dpf.post.simulation import (
    HarmonicMechanicalSimulation,
    ModalMechanicalSimulation,
    StaticMechanicalSimulation,
    TransientMechanicalSimulation,
)


def test_load_simulation_static_mechanical(simple_bar):
    simulation = dpf.load_simulation(simple_bar)
    assert type(simulation) == StaticMechanicalSimulation


def test_load_simulation_transient_mechanical(transient_rst):
    simulation = dpf.load_simulation(transient_rst)
    assert type(simulation) == TransientMechanicalSimulation


def test_load_simulation_modal_mechanical(modalallkindofcomplexity):
    simulation = dpf.load_simulation(modalallkindofcomplexity)
    assert type(simulation) == ModalMechanicalSimulation


def test_load_simulation_harmonic_mechanical(complex_model):
    simulation = dpf.load_simulation(complex_model)
    assert type(simulation) == HarmonicMechanicalSimulation
