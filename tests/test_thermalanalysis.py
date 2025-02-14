# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
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

import numpy as np
import pytest

from ansys.dpf import post
from ansys.dpf.post import errors as dpf_errors
from ansys.dpf.post.common import _PhysicsType
from ansys.dpf.post.temperature import Temperature


def test_thermal_steadystate(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    assert solution._model.metadata.result_info.physics_type == _PhysicsType.thermal
    temp = solution.temperature()
    assert isinstance(temp, Temperature)
    s = temp.scalar
    assert s.num_fields == 1
    assert np.isclose(s[0].data[23], 29.6247641917003)
    assert s[0].location == post.locations.nodal

    # with dpf.core operator
    from ansys.dpf import core

    op = core.Operator("TEMP")
    # op.inputs.requested_location.connect(core.locations.nodal)
    op.inputs.data_sources.connect(core.DataSources(rth_steady_state))
    fc = op.outputs.fields_container()
    assert len(fc) == s.num_fields
    assert fc[0].location == s[0].location
    assert len(fc[0].data) == len(s[0].data)
    assert np.allclose(s[0].data.tolist(), fc[0].data.tolist())
    comp = core.operators.logic.identical_fc()
    comp.inputs.fields_containerA.connect(fc)
    comp.inputs.fields_containerB.connect(s.result_fields_container)
    out = comp.outputs.boolean()
    assert out == True


to_return = "node scoping and element scoping returns the same"


def test_steadystate_nodscoping(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    temp = solution.temperature(node_scoping=[2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    assert np.allclose(s[0].data.tolist(), [41.8243869, 40.29943406])


@pytest.mark.skipif(True, reason="element scoping not available with thermal results.")
def test_steadystate_elemscoping(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    temp = solution.temperature(element_scoping=[2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    # assert np.allclose(s[0].data.tolist(), [41.8243869 , 40.29943406])
    raise Exception(to_return)


def test_steadystate_nodlocation(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    temp = solution.temperature(location=post.locations.nodal)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal


def test_steadystate_elemlocation(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    with pytest.raises(dpf_errors.NodalLocationError):
        temp = solution.temperature(location=post.locations.elemental)


def test_steadystate_elemnodallocation(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    with pytest.raises(dpf_errors.NodalLocationError):
        temp = solution.temperature(location=post.locations.elemental_nodal)


def test_steadystate_timescoping(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    temp = solution.temperature(time_scoping=1)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 30.387240610202973)


def test_steadystate_time(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    temp = solution.temperature(time=1.0)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 30.387240610202973)


def test_steadystate_set(rth_steady_state):
    solution = post.load_solution(rth_steady_state)
    temp = solution.temperature(set=1)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 30.387240610202973)


def test_thermal_transient(rth_transient):
    solution = post.load_solution(rth_transient)
    assert solution._model.metadata.result_info.physics_type == _PhysicsType.thermal
    temp = solution.temperature()
    assert isinstance(temp, Temperature)
    s = temp.scalar
    assert s.num_fields == 1
    assert np.isclose(s[0].data[23], 22.159122145252788)


def test_transient_nodscoping(rth_transient):
    solution = post.load_solution(rth_transient)
    temp = solution.temperature(node_scoping=[2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    assert np.allclose(s[0].data.tolist(), [27.01872925, 25.61222966])


@pytest.mark.skipif(True, reason="element scoping not available with thermal results.")
def test_transient_elemscoping(rth_transient):
    solution = post.load_solution(rth_transient)
    temp = solution.temperature(element_scoping=[2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    # assert np.allclose(s[0].data.tolist(), [27.01872925, 25.61222966])
    raise Exception(to_return)


def test_transient_nodlocation(rth_transient):
    solution = post.load_solution(rth_transient)
    temp = solution.temperature(location=post.locations.nodal)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal


def test_transient_elemlocation(rth_transient):
    solution = post.load_solution(rth_transient)
    with pytest.raises(dpf_errors.NodalLocationError):
        temp = solution.temperature(location=post.locations.elemental)


def test_transient_elemnodallocation(rth_transient):
    solution = post.load_solution(rth_transient)
    with pytest.raises(dpf_errors.NodalLocationError):
        temp = solution.temperature(location=post.locations.elemental_nodal)


def test_transient_timescoping(rth_transient):
    solution = post.load_solution(rth_transient)
    temp = solution.temperature(time_scoping=[2, 15])
    s = temp.scalar
    assert s.num_fields == 2
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 21.999999999992536)


def test_transient_time(rth_transient):
    solution = post.load_solution(rth_transient)
    temp = solution.temperature(set=2)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 21.999999999992536)


def test_transient_set(rth_transient):
    solution = post.load_solution(rth_transient)
    temp = solution.temperature(time=0.0223)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 21.999999999992323)


def test_heat_flux(rth_transient):
    solution = post.load_solution(rth_transient)
    hf = solution.heat_flux(location=post.locations.elemental)
    s = hf.vector
    assert len(s[0].data) == 784
    assert s[0].location == post.locations.elemental
    assert np.allclose(s[0].data[24], [-3.85171006e-10, -9.35413524e-10, 1.81041315e03])

    # with dpf.core operator
    from ansys.dpf import core

    op = core.Operator("TF")
    op.inputs.requested_location.connect(core.locations.elemental)
    op.inputs.data_sources.connect(core.DataSources(rth_transient))
    fc = op.outputs.fields_container()
    assert len(fc) == s.num_fields
    assert fc[0].location == s[0].location
    assert len(fc[0].data) == len(s[0].data)
    assert np.allclose(s[0].data.tolist(), fc[0].data.tolist())
    comp = core.operators.logic.identical_fc()
    comp.inputs.fields_containerA.connect(fc)
    comp.inputs.fields_containerB.connect(s.result_fields_container)
    out = comp.outputs.boolean()
    assert out == True
