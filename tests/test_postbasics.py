import os

from ansys import dpf
from ansys.dpf import post
from ansys.dpf.post.common import _AnalysisType
from ansys.dpf.post.static_analysis import StaticAnalysisSolution
from ansys.dpf.post.modal_analysis import ModalAnalysisSolution
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
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