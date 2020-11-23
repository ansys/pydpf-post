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
    
    
def test_num_fields():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    assert isinstance(disp, ResultData)
    assert disp.num_fields() == 1
    
    
def test_data_at_field():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    data = disp.data_at_field(0)
    assert data.__len__() == 15113
    assert data[0].__len__() == 3
    
    
def test_field_getitem():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    field = disp[0]
    assert isinstance(field, dpf.core.Field)
    assert field.location == locations.nodal
    
    
def test_max():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    max_val = disp.max()
    assert max_val.__len__() == 3
    assert max_val.data.__len__() == 1
    assert max_val.data[0].__len__() == 3


def test_min():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    min_val = disp.min()
    assert min_val.__len__() == 3
    assert min_val.data.__len__() == 1
    assert min_val.data[0].__len__() == 3
    

def test_maxdata():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    val = disp.max_data()
    assert val.__len__() == 1
    assert val[0].__len__() == 3
    
    
def test_mindata():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    val = disp.min_data()
    assert val.__len__() == 1
    assert val[0].__len__() == 3
    
    
def test_maxdata_at_field():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    val = disp.max_data_at_field(0)
    assert val.__len__() == 3
    
    
def test_min_data_at_field():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement()
    val = disp.min_data_at_field(0)
    assert val.__len__() == 3
    