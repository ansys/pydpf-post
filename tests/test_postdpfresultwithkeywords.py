import os
from ansys import dpf
from ansys.dpf import post
from ansys.dpf.core.common import locations
import numpy as np


if 'AWP_UNIT_TEST_FILES' in os.environ:
    unit_test_files = os.environ['AWP_UNIT_TEST_FILES']
else:
    raise KeyError('Please add the location of the DataProcessing '
                   'test files "AWP_UNIT_TEST_FILES" to your env')
    

TEST_FILE_PATH_RST = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'allKindOfComplexity.rst')


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
def test_displacement_with_scoping():
    result = post.result(TEST_FILE_PATH_RST)
    #scoping as array
    disp = result.nodal_displacement(node_scoping = [1, 2])
    data = disp.data_at_field(0)
    assert data.__len__() == 2
    assert data[0].__len__() == 3
    #scoping as dpf.Scoping
    scop = dpf.core.Scoping()
    scop.ids = [1, 2]
    scop.location = locations.nodal
    disp2 = result.nodal_displacement(node_scoping = scop)
    data2 = disp2.data_at_field(0)
    assert data2.__len__() == 2
    assert data2[0].__len__() == 3
    #values comparison
    assert np.allclose(data, data2)
            

def test_node_stress_with_scoping():
    result = post.result(TEST_FILE_PATH_RST)
    #scoping as array
    disp = result.nodal_stress(element_scoping = [1, 34])
    data = disp.data_at_field(0)
    assert data.__len__() == 40
    assert data[0].__len__() == 6
    #scoping as dpf.Scoping
    scop = dpf.core.Scoping()
    scop.location = locations.elemental
    scop.ids = [1, 34]
    disp2 = result.nodal_stress(element_scoping = scop)
    data2 = disp2.data_at_field(0)
    assert data2.__len__() == 40
    assert data2[0].__len__() == 6
    #values comparison
    assert np.allclose(data, data2)
    
    
def test_elem_stress_with_scoping():
    result = post.result(TEST_FILE_PATH_RST)
    #scoping as array
    disp = result.elemental_stress(element_scoping = [1, 34])
    data = disp.data_at_field(0)
    assert data.__len__() == 16
    assert data[0].__len__() == 6
    #scoping as dpf.Scoping
    scop = dpf.core.Scoping()
    scop.location = locations.elemental
    scop.ids = [1, 34]
    disp2 = result.elemental_stress(element_scoping = scop)
    data2 = disp2.data_at_field(0)
    assert data2.__len__() == 16
    assert data2[0].__len__() == 6
    #values comparison
    assert np.allclose(data, data2)
    
    
def test_disp_with_component_subresult():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement(subresult="Y")
    assert disp._result_operator.name == "UY"
    assert disp.num_fields() == 1
    data = disp.data_at_field(0)
    assert data.__len__() == 15113
    assert data[0] == 5.130250313479703e-06


def test_stress_with_component_subresult():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.elemental_stress(subresult="YZ")
    assert stress._result_operator.name == "SYZ"
    assert stress.num_fields() == 1
    data = stress.data_at_field(0)
    assert data.__len__() == 40016
    assert data[1] == 1.0216815465593042e-10


def test_stress_with_invariant_subresult():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.elemental_stress(subresult="3")
    assert stress._result_operator.name == "S3"
    assert stress.num_fields() == 2
    data = stress.data_at_field(0)
    assert data.__len__() == 720
    assert data[1] == -4721842.179373354


def test_von_mises_stress():
    raise Exception("Not impl. yet")
            

    
    