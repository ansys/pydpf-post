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

from ansys.dpf.core import examples
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
class TestPhase:
    @fixture
    def fluid_simulation(self):
        fluid_example_files = examples.download_fluent_axial_comp()
        ds = dpf.DataSources()
        ds.set_result_file_path(
            fluid_example_files["cas"][0],
            key="cas",
        )
        ds.add_file_path(
            fluid_example_files["dat"][0],
            key="dat",
        )
        return post.FluidSimulation(ds)  # noqa

    def test_phase(self):
        phase = post.Phase(name="phase_test", id=1)
        assert isinstance(phase, post.Phase)
        assert phase.name == "phase_test"
        assert phase.id == 1
        assert repr(phase) == "Phase<name: 'phase_test', id=1>"

    def test_phases(self, fluid_simulation):
        phases = post.PhasesDict(simulation=fluid_simulation)
        assert isinstance(phases, post.PhasesDict)
        assert len(phases) == 1
        assert isinstance(phases[1], post.Phase)
        assert repr(phases) == "{Phase<name: 'phase-1', id=1>, }"
        assert str(phases) == "1 phases available\n{Phase<name: 'phase-1', id=1>, }"
        with pytest.raises(ValueError):
            _ = phases["toto"]
        with pytest.raises(ValueError):
            _ = phases[32]
        with pytest.raises(ValueError):
            _ = phases[24.6]
        for phase in phases:
            assert isinstance(phase, post.Phase)
        for phase in phases:
            assert isinstance(phase, post.Phase)
        _ = list(phases)
