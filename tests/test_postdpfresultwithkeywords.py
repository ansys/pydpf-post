import os
from ansys import dpf
from ansys.dpf import post
from ansys.dpf.core.common import locations
import numpy as np
import unittest


if 'AWP_UNIT_TEST_FILES' in os.environ:
    unit_test_files = os.environ['AWP_UNIT_TEST_FILES']
else:
    raise KeyError('Please add the location of the DataProcessing '
                   'test files "AWP_UNIT_TEST_FILES" to your env')
    

TEST_FILE_PATH_RST = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'allKindOfComplexity.rst')

TRANSIENT_FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'expansion', 
                                   'msup', 'Transient', 'plate1','file.rst')


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
    
    
def test_elem_nodal_stress_with_scoping():
    result = post.result(TEST_FILE_PATH_RST)
    #scoping as array
    disp = result.elemental_nodal_stress(element_scoping = [1, 34])
    data = disp.data_at_field(0)
    assert data.__len__() == 16
    assert data[0].__len__() == 6
    #scoping as dpf.Scoping
    scop = dpf.core.Scoping()
    scop.location = locations.elemental
    scop.ids = [1, 34]
    disp2 = result.elemental_nodal_stress(element_scoping = scop)
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
    stress = result.elemental_nodal_stress(subresult="YZ")
    assert stress._result_operator.name == "SYZ"
    assert stress.num_fields() == 1
    data = stress.data_at_field(0)
    assert data.__len__() == 40016
    assert data[1] == 1.0216815465593042e-10


def test_stress_with_invariant_subresult():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.elemental_nodal_stress(subresult="3")
    assert stress._result_operator.name == "S3"
    assert stress.num_fields() == 2
    data = stress.data_at_field(0)
    assert data.__len__() == 720
    assert data[1] == -4721842.179373354


def test_von_mises_stress():
    raise Exception("Not impl. yet")
    
    
def test_groupingelshape_nodallocation():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement(grouping = post.grouping.by_el_shape)
    assert disp.num_fields() == 4
    assert disp.result_fields_container.get_label_space(3) == {'elshape': 3, 'time': 1}
    assert disp.data_at_field(0).__len__() == 93950
    assert disp.data_at_field(1).__len__() == 6748
    assert disp.data_at_field(2).__len__() == 36
    assert disp.data_at_field(3).__len__() == 4
    assert disp.data_at_field(2)[0][0] == 5.523488975819807e-20
    assert disp[0].location == locations.nodal


def test_groupingelshape_elemlocation():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.elemental_stress(grouping = post.grouping.by_el_shape)
    assert stress.num_fields() == 4
    assert stress.result_fields_container.get_label_space(3) == {'elshape': 3, 'time': 1}
    assert stress.data_at_field(0).__len__() == 609
    assert stress.data_at_field(1).__len__() == 27156
    assert stress.data_at_field(2).__len__() == 0
    assert stress.data_at_field(3).__len__() == 0
    assert stress.data_at_field(1)[0][0] == -2984789.100540372
    assert stress[0].location == locations.elemental
    

def test_groupingmat_nodallocation():
    result = post.result(TEST_FILE_PATH_RST)
    disp = result.nodal_displacement(grouping = post.grouping.by_material)
    assert disp.num_fields() == 11
    assert disp[0].__len__() == 23016
    assert disp[2].__len__() == 1848
    assert disp.data_at_field(2)[0][2] == -6.649053654123576e-07
    assert disp.result_fields_container.get_label_space(3) == {'time': 1, 'mat': 10}
    

def test_groupingmat_elemlocation():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.elemental_stress(grouping = post.grouping.by_material)
    assert stress.num_fields() == 11
    # assert stress[0].__len__() == 0
    # assert stress[5].__len__() == 156762
    # assert stress.data_at_field(5)[0][2] == 2255995657.8043337
    assert stress.result_fields_container.get_label_space(3) == {'time': 1, 'mat': 4}
    assert stress[0].location == locations.elemental
    

def test_mapdlgrouping_nodallocation():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.nodal_displacement(mapdl_grouping = 186)
    try:
        stress.num_fields()
    except:
        assert True
    

def test_maplgrouping_elemlocation():
    result = post.result(TEST_FILE_PATH_RST)
    stress = result.elemental_stress(mapdl_grouping = 186)
    assert stress.num_fields() == 1
    assert stress.result_fields_container.get_label_space(0) == {'time': 1}
    assert stress.data_at_field(0).__len__() == 343
    assert stress.data_at_field(0)[41][2] == -323198.184976747
    assert stress[0].location == locations.elemental
    

def test_set_keyword():
    result = post.result(TRANSIENT_FILE_PATH)
    disp = result.nodal_displacement(set = 3)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 3}
    assert disp.data_at_field(0)[2][2] == 2.3955190605044603e-05
    
    
class TestCase(unittest.TestCase):
    def test_both_set_time(self):
        result = post.result(TRANSIENT_FILE_PATH)
        self.assertRaises(Exception, result.nodal_displacement, set=3, time=0.01)
        try:
            result.nodal_displacement(set = 3, time = 0.01)
        except Exception as e:
            message = "'time' and 'set' keyword can not be used simultaneously."
            e2 = Exception(message)
            assert e.args == e2.args
            assert type(e) == type(e2)
            

def test_time_keyword_in_frequencies():
    result = post.result(TRANSIENT_FILE_PATH)
    disp = result.nodal_displacement(time=0.06)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 6}
    assert disp.data_at_field(0)[2][2] == 6.449354759605568e-05
    disp = result.nodal_displacement(time=0.02)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 2}
    assert disp.data_at_field(0)[40][2] == -9.555678764252377e-06
    disp = result.nodal_displacement(time=0.14)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 14}
    assert disp.data_at_field(0)[22][2] == -5.9753488295405e-06
    disp = result.nodal_displacement(time=0.15)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 15}
    assert disp.data_at_field(0)[101][2] == 1.2454347438346573e-05
    disp = result.nodal_displacement(time=0.2)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 20}
    assert disp.data_at_field(0)[345][2] == 6.931130871751968e-05


def test_time_keyword_not_in_frequencies():
    result = post.result(TRANSIENT_FILE_PATH)
    disp = result.nodal_displacement(time=0.061)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 0}
    assert disp.data_at_field(0)[2][2] == 6.466312449668174e-05
    disp = result.nodal_displacement(time=0.023)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 0}
    assert disp.data_at_field(0)[40][2] == -1.3341949773184135e-05
    disp = result.nodal_displacement(time=0.1499)
    assert disp.num_fields() == 1
    assert disp.result_fields_container.get_label_space(0) == {'time': 0}
    assert disp.data_at_field(0)[22][2] == -1.7795334817918245e-05
    
    
def test_time_scoping_keyword():
    result = post.result(TRANSIENT_FILE_PATH)
    disp = result.nodal_displacement()
    assert disp.num_fields() == 1
    disp1 = result.nodal_displacement(time_scoping=[1,2,4])
    assert disp1.num_fields() == 3
    assert disp1.result_fields_container.get_label_space(0) == {'time': 1}
    assert disp1.data_at_field(0)[40][2] == -2.0115581116044217e-06
    # disp2 = result.nodal_displacement(time_scoping=np.array([1, 2, 4]))
    # assert disp2.num_fields() == 3
    # assert disp2.result_fields_container.get_label_space(0) == {'time': 1}
    # assert disp2.data_at_field(0)[40][2] == -2.0115581116044217e-06
    scop = dpf.core.Scoping()
    scop.ids = [1, 2, 4]
    disp3 = result.nodal_displacement(time_scoping=scop)
    assert disp3.num_fields() == 3
    assert disp3.result_fields_container.get_label_space(0) == {'time': 1}
    assert disp3.data_at_field(0)[40][2] == -2.0115581116044217e-06
    
    
    