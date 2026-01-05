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

import numpy as np

from ansys.dpf import post
from ansys.dpf.post.displacement import ComplexDisplacement
from ansys.dpf.post.dpf_solution import DpfMecanicComplexSolution
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.result_data import ResultData


def test_displacement_amplitude_verbose_api(complex_model):
    result = post.load_solution(complex_model)
    assert isinstance(result, HarmonicAnalysisSolution)
    assert isinstance(result, DpfMecanicComplexSolution)
    ampl = result.misc.nodal_displacement_amplitude()
    assert isinstance(ampl, ResultData)
    assert ampl.num_fields == 1
    l = ampl.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert np.isclose(l[2][1], 2.863712949981395e-10)


def test_displacement_amplitude(complex_model):
    sol = post.load_solution(complex_model)
    assert isinstance(sol, HarmonicAnalysisSolution)
    assert isinstance(sol, DpfMecanicComplexSolution)
    complex_disp = sol.displacement()
    assert isinstance(complex_disp, ComplexDisplacement)
    ampl = complex_disp.vector_amplitude
    assert isinstance(ampl, ResultData)
    assert ampl.num_fields == 1
    l = ampl.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert np.isclose(l[2][1], 2.863712949981395e-10)


def test_displacement_at_phase_verbose_api(complex_model):
    result = post.load_solution(complex_model)
    assert isinstance(result, HarmonicAnalysisSolution)
    assert isinstance(result, DpfMecanicComplexSolution)
    disp_at_phase = result.misc.nodal_displacement(phase=41.0)
    assert isinstance(disp_at_phase, ResultData)
    assert disp_at_phase.num_fields == 1
    l = disp_at_phase.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert np.isclose(l[2][1], -2.1606921923925902e-10)


def test_displacement_at_phase(complex_model):
    result = post.load_solution(complex_model)
    assert isinstance(result, HarmonicAnalysisSolution)
    assert isinstance(result, DpfMecanicComplexSolution)
    complex_disp = result.displacement()
    assert isinstance(complex_disp, ComplexDisplacement)
    disp_at_phase = complex_disp.vector_at_phase(phase=41.0)
    assert isinstance(disp_at_phase, ResultData)
    assert disp_at_phase.num_fields == 1
    l = disp_at_phase.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert np.isclose(l[2][1], -2.1606921923925902e-10)


def test_has_complex_result(complex_model):
    result = post.load_solution(complex_model)
    assert result.has_complex_result()


# def test_is_complex_result_verbose_api(complex_model):
#     result = post.load_solution(complex_model)
#     disp = result.misc.nodal_displacement()
#     assert disp.num_fields == 2


# def test_is_complex_result_verbose_api_modal(modalallkindofcomplexity):
#     result = post.load_solution(modalallkindofcomplexity)
#     disp = result.misc.nodal_displacement()
#     assert disp.num_fields == 1


def test_is_complex_result(complex_model):
    result = post.load_solution(complex_model)
    compl_disp = result.displacement()
    assert compl_disp.has_complex_frequencies()
    disp = compl_disp.vector
    assert disp.num_fields == 2
