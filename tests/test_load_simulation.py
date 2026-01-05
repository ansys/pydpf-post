# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import ansys.dpf.core as dpf
import pytest

from ansys.dpf import post
from ansys.dpf.post.common import AvailableSimulationTypes
from ansys.dpf.post.fluid_simulation import FluidSimulation
from ansys.dpf.post.harmonic_mechanical_simulation import HarmonicMechanicalSimulation
from ansys.dpf.post.modal_mechanical_simulation import ModalMechanicalSimulation
from ansys.dpf.post.static_mechanical_simulation import StaticMechanicalSimulation
from ansys.dpf.post.transient_mechanical_simulation import TransientMechanicalSimulation
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
)


def test_load_simulation_static_mechanical(simple_bar, complex_model):
    simulation = post.load_simulation(data_sources=simple_bar)
    assert isinstance(simulation, StaticMechanicalSimulation)
    simulation = post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.static_mechanical,
    )
    assert isinstance(simulation, StaticMechanicalSimulation)


def test_load_simulation_transient_mechanical(plate_msup, complex_model):
    simulation = post.load_simulation(data_sources=plate_msup)
    assert isinstance(simulation, TransientMechanicalSimulation)
    simulation = post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )
    assert isinstance(simulation, TransientMechanicalSimulation)


def test_load_simulation_modal_mechanical(modalallkindofcomplexity, complex_model):
    simulation = post.load_simulation(data_sources=modalallkindofcomplexity)
    assert isinstance(simulation, ModalMechanicalSimulation)
    simulation = post.load_simulation(
        data_sources=complex_model,
        simulation_type=AvailableSimulationTypes.transient_mechanical,
    )
    assert isinstance(simulation, TransientMechanicalSimulation)


def test_load_simulation_harmonic_mechanical(complex_model, simple_bar):
    simulation = post.load_simulation(data_sources=complex_model)
    assert isinstance(simulation, HarmonicMechanicalSimulation)
    simulation = post.load_simulation(
        data_sources=simple_bar,
        simulation_type=AvailableSimulationTypes.harmonic_mechanical,
    )
    assert isinstance(simulation, HarmonicMechanicalSimulation)


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
def test_load_simulation_static_fluid(fluid_fluent_elbow_steady_state, simple_bar):
    ds = dpf.DataSources()
    ds.set_result_file_path(fluid_fluent_elbow_steady_state["cas"][0], key="cas")
    ds.add_file_path(
        fluid_fluent_elbow_steady_state["dat"][0],
        key="dat",
    )
    simulation = post.load_simulation(ds)
    assert isinstance(simulation, FluidSimulation)
    simulation = post.load_simulation(
        data_sources=simple_bar,
        simulation_type=AvailableSimulationTypes.steady_fluid,
    )
    assert isinstance(simulation, FluidSimulation)


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
def test_load_simulation_transient_fluid(fluid_fluent_elbow_transient, simple_bar):
    ds = dpf.DataSources()
    ds.set_result_file_path(fluid_fluent_elbow_transient["cas"][0], key="cas")
    ds.add_file_path(
        fluid_fluent_elbow_transient["dat"][0],
        key="dat",
    )
    simulation = post.load_simulation(ds)
    assert isinstance(simulation, FluidSimulation)
    simulation = post.load_simulation(
        data_sources=simple_bar,
        simulation_type=AvailableSimulationTypes.unsteady_fluid,
    )
    assert isinstance(simulation, FluidSimulation)


def test_load_simulation_raise_simulation_type(simple_bar):
    with pytest.raises(ValueError, match="is not a recognized simulation type"):
        _ = post.load_simulation(
            data_sources=simple_bar,
            simulation_type="test",
        )


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0,
    reason="Available starting DPF 4.0",
)
def test_load_simulation_with_server(simple_bar, grpc_server):
    simulation = post.load_simulation(data_sources=simple_bar, server=grpc_server)
    assert isinstance(simulation, StaticMechanicalSimulation)
    assert simulation._model._server != dpf.SERVER
    assert simulation._model._server == grpc_server
