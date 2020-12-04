import os
from ansys import dpf
from ansys.dpf import post
from ansys.dpf.post.result_data import ResultData
from ansys.dpf.core.common import locations


if 'AWP_UNIT_TEST_FILES' in os.environ:
    unit_test_files = os.environ['AWP_UNIT_TEST_FILES']
else:
    raise KeyError('Please add the location of the DataProcessing '
                   'test files "AWP_UNIT_TEST_FILES" to your env')
    

TEST_FILE_PATH_RST = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'allKindOfComplexity.rst')


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
def test_get_result_info():
    result = post.solution(TEST_FILE_PATH_RST)
    res_info = result.get_result_info()
    assert res_info.analysis_type
    assert res_info.unit_system
    assert res_info.physics_type
    assert res_info.available_results
    

def test_get_result_data_function_of_operator_no_keyword():
    result = post.solution(TEST_FILE_PATH_RST)
    result_data = result._get_result_data_function_of_operator("U", result, result._data_sources)
    assert isinstance(result_data, ResultData)
    
    
def test_get_result_data_function_of_operator_ns():
    result = post.solution(TEST_FILE_PATH_RST)
    result_data = result._get_result_data_function_of_operator("U", result, result._data_sources, named_selection="SELECTION")
    assert isinstance(result_data, ResultData)


def test_get_result_data_function_of_operator_location():
    result = post.solution(TEST_FILE_PATH_RST)
    result_data = result._get_result_data_function_of_operator("U", result, result._data_sources, location="Elemental")
    assert isinstance(result_data, ResultData)


def test_get_result_data_function_of_operator_node_scop():
    result = post.solution(TEST_FILE_PATH_RST)
    result_data = result._get_result_data_function_of_operator("U", result, result._data_sources, node_scoping=[1, 2, 3])
    assert isinstance(result_data, ResultData)


def test_get_result_data_function_of_operator_elem_scop():
    result = post.solution(TEST_FILE_PATH_RST)
    result_data = result._get_result_data_function_of_operator("U", result, result._data_sources, element_scoping=[1, 2, 3])
    assert isinstance(result_data, ResultData)


def test_get_result_data_function_of_operator_bothscop():
    result = post.solution(TEST_FILE_PATH_RST)
    try:
        result._get_result_data_function_of_operator("U", result, result._data_sources, node_scoping=[1, 2, 3], element_scoping=[1, 2, 3])
    except Exception as e:
        message = "Impossible to use both element_scoping and node_scoping."
        e2 = Exception(message)
        assert e.args == e2.args
        assert type(e) == type(e2)


def test_get_result_data_function_of_operator_nophase():
    result = post.solution(TEST_FILE_PATH_RST)
    try:
        result._get_result_data_function_of_operator("U", result, result._data_sources, phase=30.)
    except Exception as e:
        message = "Phase key-word argument can be used when the analysis types implies complex result (Harmonic analysis, Modal analysis...)."
        e2 = Exception(message)
        assert e.args == e2.args
        assert type(e) == type(e2)


def test_check_elemental_location():
    result = post.solution(TEST_FILE_PATH_RST)
    try:
        result.nodal_displacement(location="Elemental")
    except Exception as e:
        message = "Only a Nodal location can be used with a nodal result."
        e2 = Exception(message)
        assert e.args == e2.args
        assert type(e) == type(e2)
        

def test_check_nodal_location():
    result = post.solution(TEST_FILE_PATH_RST)
    try:
        result.elemental_stress(location="nodal")
    except Exception as e:
        message = "Only an Elemental location can be used with an elemental result."
        e2 = Exception(message)
        assert e.args == e2.args
        assert type(e) == type(e2)


def test_nodal_displacement():
    result = post.solution(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    assert isinstance(disp, ResultData)
    assert disp.num_fields() == 1
    data = disp.data_at_field(0)
    assert len(data) == 15113
    assert len(data[0]) == 3
    field = disp[0]
    assert field.location == locations.nodal
    
    
def test_nodal_stress():
    result = post.solution(TEST_FILE_PATH_RST)
    stress = result.nodal_stress()
    assert isinstance(stress, ResultData)
    assert stress.num_fields() == 2
    data = stress.data_at_field(0)
    assert len(data) == 720
    assert len(data[0]) == 6
    data = stress.data_at_field(1)
    assert len(data) == 14826
    assert len(data[0]) == 6
    field = stress[0]
    assert field.location == locations.nodal


def test_elemental_stress():
    result = post.solution(TEST_FILE_PATH_RST)
    stress = result.elemental_stress()
    assert isinstance(stress, ResultData)
    assert stress.num_fields() == 2
    data = stress.data_at_field(0)
    assert len(data) == 609
    assert len(data[0]) == 6
    data = stress.data_at_field(1)
    assert len(data) == 9052
    assert len(data[0]) == 6
    field = stress[0]
    assert field.location == locations.elemental

    
def test_elemental_nodal_stress():
    result = post.solution(TEST_FILE_PATH_RST)
    stress = result.elemental_nodal_stress()
    assert isinstance(stress, ResultData)
    assert stress.num_fields() == 1
    data = stress.data_at_field(0)
    assert len(data) == 40016
    assert len(data[0]) == 6
    field = stress[0]
    assert field.location == locations.elemental_nodal
    
        
        