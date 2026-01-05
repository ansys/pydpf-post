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

from ansys.dpf import post
from ansys.dpf.post.common import _AnalysisType
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.modal_analysis import ModalAnalysisSolution
from ansys.dpf.post.static_analysis import StaticAnalysisSolution
from ansys.dpf.post.transient_analysis import TransientAnalysisSolution


def test_call_result_object_static(simple_bar):
    result = post.load_solution(simple_bar)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.static
    assert isinstance(result, StaticAnalysisSolution)


def test_call_result_object_modal(modalallkindofcomplexity):
    result = post.load_solution(modalallkindofcomplexity)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.modal
    assert isinstance(result, ModalAnalysisSolution)


def test_call_result_object_harmonic(complex_model):
    result = post.load_solution(complex_model)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.harmonic
    assert isinstance(result, HarmonicAnalysisSolution)


def test_call_result_object_transient(plate_msup):
    result = post.load_solution(plate_msup)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.transient
    assert isinstance(result, TransientAnalysisSolution)
