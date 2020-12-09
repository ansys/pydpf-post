import os
from ansys import dpf
from ansys.dpf import post
from ansys.dpf.post.dpf_solution import DpfComplexSolution
from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.displacement import ComplexDisplacement


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
    
    
def test_displacement_amplitude_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    assert isinstance(result, HarmonicAnalysisSolution)
    assert isinstance(result, DpfComplexSolution)
    ampl = result.misc.nodal_displacement_amplitude()
    assert isinstance(ampl, ResultData)
    assert ampl.num_fields == 1
    l = ampl.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert l[2][1] == 2.863712949981395e-10
    

def test_displacement_amplitude():
    sol = post.load_solution(TEST_FILE_PATH_RST)
    assert isinstance(sol, HarmonicAnalysisSolution)
    assert isinstance(sol, DpfComplexSolution)
    complex_disp = sol.displacement()
    assert isinstance(complex_disp, ComplexDisplacement)
    ampl = complex_disp.vector_amplitude
    assert isinstance(ampl, ResultData)
    assert ampl.num_fields == 1
    l = ampl.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert l[2][1] == 2.863712949981395e-10
    
    
def test_displacement_at_phase_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    assert isinstance(result, HarmonicAnalysisSolution)
    assert isinstance(result, DpfComplexSolution)
    disp_at_phase = result.misc.nodal_displacement(phase = 41.)
    assert isinstance(disp_at_phase, ResultData)
    assert disp_at_phase.num_fields == 1
    l = disp_at_phase.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert l[2][1] == -2.1606921923925902e-10
    

def test_displacement_at_phase():
    result = post.load_solution(TEST_FILE_PATH_RST)
    assert isinstance(result, HarmonicAnalysisSolution)
    assert isinstance(result, DpfComplexSolution)
    complex_disp = result.displacement()
    assert isinstance(complex_disp, ComplexDisplacement)
    disp_at_phase = complex_disp.vector_at_phase(phase = 41.)
    assert isinstance(disp_at_phase, ResultData)
    assert disp_at_phase.num_fields == 1
    l = disp_at_phase.get_data_at_field(0)
    assert len(l) == 4802
    assert len(l[2]) == 3
    assert l[2][1] == -2.1606921923925902e-10
    
    
def test_has_complex_result():
    result = post.load_solution(TEST_FILE_PATH_RST)
    assert result.has_complex_result()
    

def test_is_complex_result_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    assert disp.num_fields == 2
    assert disp.is_complex_result()
    result = post.load_solution(MODAL_FILE_PATH)
    disp = result.misc.nodal_displacement()
    assert disp.num_fields == 1
    assert disp.is_complex_result()
    
    
def test_is_complex_result():
    result = post.load_solution(TEST_FILE_PATH_RST)
    compl_disp = result.displacement()
    disp = compl_disp.vector
    assert disp.num_fields == 2
    assert disp.is_complex_result()
    result = post.load_solution(MODAL_FILE_PATH)
    disp = result.misc.nodal_displacement()
    assert disp.num_fields == 1
    assert disp.is_complex_result()