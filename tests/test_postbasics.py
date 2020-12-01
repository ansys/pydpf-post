import os

from ansys import dpf
from ansys.dpf import post
from ansys.dpf.post.common import _AnalysisType
from ansys.dpf.post.static_analysis import StaticAnalysisResult
from ansys.dpf.post.modal_analysis import ModalAnalysisResult
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisResult
from ansys.dpf.post.transient_analysis import TransientAnalysisResult


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
if 'AWP_UNIT_TEST_FILES' in os.environ:
    unit_test_files = os.environ['AWP_UNIT_TEST_FILES']
else:
    raise KeyError('Please add the location of the DataProcessing '
                   'test files "AWP_UNIT_TEST_FILES" to your env')
    
    
STATIC_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'ASimpleBar.rst')

MODAL_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'modal_allKindOfComplexity.rst')

HARMONIC_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'rth_operators',
                              'fileComplex.rst') #with complex results

TRANSIENT_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'expansion', 
                                   'msup', 'Transient', 'plate1','file.rst')


def test_build_docs():
    post.build_docs()
    assert os.path.exists("dataProcessingDoc.html")
    os.remove("dataProcessingDoc.html")
    assert not os.path.exists("dataProcessingDoc.html")
    
    
def test_call_result_object_static():
    result = post.result(STATIC_FILE_PATH)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.static
    assert isinstance(result, StaticAnalysisResult)
    
    
def test_call_result_object_modal():
    result = post.result(MODAL_FILE_PATH)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.modal
    assert isinstance(result, ModalAnalysisResult)
    
    
def test_call_result_object_harmonic():
    result = post.result(HARMONIC_FILE_PATH)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.harmonic
    assert isinstance(result, HarmonicAnalysisResult)
    
    
def test_call_result_object_transient():
    result = post.result(TRANSIENT_FILE_PATH)
    assert result._model.metadata.result_info.analysis_type == _AnalysisType.transient
    assert isinstance(result, TransientAnalysisResult)