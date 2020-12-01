import os
import unittest
from ansys import dpf
from ansys.dpf import post
from ansys.dpf.post.dpf_result import DpfComplexResult
from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisResult


if 'AWP_UNIT_TEST_FILES' in os.environ:
    unit_test_files = os.environ['AWP_UNIT_TEST_FILES']
else:
    raise KeyError('Please add the location of the DataProcessing '
                   'test files "AWP_UNIT_TEST_FILES" to your env')
    

TEST_FILE_PATH_RST = os.path.join(unit_test_files, 'DataProcessing', 'rth_operators',
                              'fileComplex.rst')


MODAL_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'modal_allKindOfComplexity.rst')


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
def test_displacement_amplitude():
    result = post.solution(TEST_FILE_PATH_RST)
    assert isinstance(result, HarmonicAnalysisResult)
    assert isinstance(result, DpfComplexResult)
    ampl = result.nodal_displacement_amplitude()
    assert isinstance(ampl, ResultData)
    assert ampl.num_fields() == 1
    l = ampl.data_at_field(0)
    assert l.__len__() == 4802
    assert l[2].__len__() == 3
    
    
def test_displacement_at_phase():
    result = post.solution(TEST_FILE_PATH_RST)
    assert isinstance(result, HarmonicAnalysisResult)
    assert isinstance(result, DpfComplexResult)
    disp_at_phase = result.nodal_displacement(phase = 41.)
    assert isinstance(disp_at_phase, ResultData)
    assert disp_at_phase.num_fields() == 1
    l = disp_at_phase.data_at_field(0)
    assert l.__len__() == 4802
    assert l[2].__len__() == 3
    
    
def test_has_complex_result():
    result = post.solution(TEST_FILE_PATH_RST)
    assert result.has_complex_result()
    

def test_is_complex_result():
    result = post.solution(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    assert disp.num_fields() == 2
    assert disp.is_complex_result()
    result = post.solution(MODAL_FILE_PATH)
    disp = result.nodal_displacement()
    assert disp.num_fields() == 1
    assert disp.is_complex_result()
    