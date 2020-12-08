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


MODAL_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'modal_allKindOfComplexity.rst')


TRANSIENT_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'expansion', 
                                   'msup', 'Transient', 'plate1','file.rst')


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
def test_num_fields_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    assert isinstance(disp, ResultData)
    assert disp.num_fields == 1
    

def test_num_fields():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    assert isinstance(disp, ResultData)
    assert disp.num_fields == 1
    
    
def test_data_at_field_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    data = disp.get_data_at_field(0)
    assert len(data) == 15113
    assert len(data[0]) == 3
    
    
def test_data_at_field():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    data = disp.get_data_at_field(0)
    assert len(data) == 15113
    assert len(data[0]) == 3
    
    
def test_field_getitem_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    field = disp[0]
    assert isinstance(field, dpf.core.Field)
    assert field.location == locations.nodal
    
    
def test_field_getitem():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    field = disp[0]
    assert isinstance(field, dpf.core.Field)
    assert field.location == locations.nodal
    
    
def test_max_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    max_val = disp.max
    assert len(max_val) == 3
    assert len(max_val.data) == 1
    assert len(max_val.data[0]) == 3
    
    
def test_max():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    max_val = disp.max
    assert len(max_val) == 3
    assert len(max_val.data) == 1
    assert len(max_val.data[0]) == 3


def test_min_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    min_val = disp.min
    assert len(min_val) == 3
    assert len(min_val.data) == 1
    assert len(min_val.data[0]) == 3
    
    
def test_min():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    min_val = disp.min
    assert len(min_val) == 3
    assert len(min_val.data) == 1
    assert len(min_val.data[0]) == 3
    

def test_maxdata_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    val = disp.max_data
    assert len(val) == 1
    assert len(val[0]) == 3
    
    
def test_maxdata():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    val = disp.max_data
    assert len(val) == 1
    assert len(val[0]) == 3
    
    
def test_mindata_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    val = disp.min_data
    assert len(val) == 1
    assert len(val[0]) == 3
    
    
def test_mindata():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    val = disp.min_data
    assert len(val) == 1
    assert len(val[0]) == 3
    
    
def test_maxdata_at_field_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    val = disp.get_max_data_at_field(0)
    assert len(val) == 3
    
    
def test_maxdata_at_field():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    val = disp.get_max_data_at_field(0)
    assert len(val) == 3
    
    
def test_min_data_at_field_verbose_api():
    result = post.load_solution(TEST_FILE_PATH_RST)
    disp = result.misc.nodal_displacement()
    val = disp.get_min_data_at_field(0)
    assert len(val) == 3
    
    
def test_min_data_at_field():
    result = post.load_solution(TEST_FILE_PATH_RST)
    d = result.displacement()
    disp = d.vector
    val = disp.get_min_data_at_field(0)
    assert len(val) == 3
    
    
def test_get_all_labels_verbose_api():
    result = post.load_solution(MODAL_FILE_PATH)
    stress = result.misc.elemental_stress()
    l = [{'elshape': 1, 'time': 1},{'elshape': 0, 'time': 1}]
    l_comp = stress.get_all_label_spaces()
    assert l == l_comp
    
    
def test_get_all_labels():
    result = post.load_solution(MODAL_FILE_PATH)
    s = result.stress(location=post.locations.elemental)
    stress = s.tensor
    l = [{'elshape': 1, 'time': 1},{'elshape': 0, 'time': 1}]
    l_comp = stress.get_all_label_spaces()
    assert l == l_comp
    
    
def test_get_scoping_at_field_verbose_api():
    result = post.load_solution(TRANSIENT_FILE_PATH)
    disp = result.misc.nodal_displacement(time_scoping=[1, 2, 4])
    assert disp.num_fields == 3
    scop = disp.get_scoping_at_field(2)
    assert len(scop) == 393
    assert scop[2] == 95
    
    
def test_get_scoping_at_field():
    result = post.load_solution(TRANSIENT_FILE_PATH)
    d = result.displacement(time_scoping=[1, 2, 4])
    disp = d.vector
    assert disp.num_fields == 3
    scop = disp.get_scoping_at_field(2)
    assert len(scop) == 393
    assert scop[2] == 95
    
    
def test_print():
    result = post.load_solution(TRANSIENT_FILE_PATH)
    print(result)
    d = result.displacement(time_scoping=[1, 2, 4])
    print(d)
    disp = d.vector
    print(disp)
    
    