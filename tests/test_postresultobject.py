from ansys import dpf
import os
import unittest
from ansys.dpf.core import locations
from ansys.dpf import post
from ansys.dpf.post.result_data import ResultData


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
if 'AWP_UNIT_TEST_FILES' in os.environ:
    unit_test_files = os.environ['AWP_UNIT_TEST_FILES']
else:
    raise KeyError('Please add the location of the DataProcessing '
                   'test files "AWP_UNIT_TEST_FILES" to your env')
    

FILE_PATH = os.path.join(unit_test_files, 'DataProcessing', 'rst_operators',
                              'allKindOfComplexity.rst')


def test_scalar():
    solution = post.load_solution(FILE_PATH)
    scalar = solution.temperature()
    txt = scalar.__str__()
    assert txt == "Scalar object. \n\nObject properties are: \n - location = Nodal\n\nThis is a temperature object."
    value = scalar.scalar
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert value[0].data[0] == 22.


def test_vector():
    solution = post.load_solution(FILE_PATH)
    vector = solution.displacement()
    txt = vector.__str__()
    assert txt == "Vector object. \n\nObject properties are: \n - location = Nodal\n\nThis is a displacement object."
    value = vector.vector
    assert isinstance(value, ResultData)
    assert value.num_fields == 1
    assert len(value[0].data) == 15113
    assert len(value[0].data[3]) == 3
    x = vector.x
    assert isinstance(x, ResultData)
    assert x.num_fields == 1
    assert len(x[0].data) == 15113
    y = vector.y
    assert isinstance(y, ResultData)
    assert y.num_fields == 1
    assert len(y[0].data) == 15113
    z = vector.z
    assert isinstance(z, ResultData)
    assert z.num_fields == 1
    assert len(z[0].data) == 15113


def test_tensor():
    solution = post.load_solution(FILE_PATH)
    tensor = solution.stress()
    txt = tensor.__str__()
    assert txt == "Stress. \nTensor object. \n\nObject properties are: \n - location = Nodal\n"
    value = tensor.tensor
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert len(value[0].data) == 720
    assert len(value[0].data[3]) == 6
    xx = tensor.xx
    assert isinstance(xx, ResultData)
    assert xx.num_fields == 2
    assert len(xx[0].data) == 720
    yy = tensor.yy
    assert isinstance(yy, ResultData)
    assert yy.num_fields == 2
    assert len(yy[0].data) == 720
    zz = tensor.zz
    assert isinstance(zz, ResultData)
    assert zz.num_fields == 2
    assert len(zz[0].data) == 720
    xy = tensor.xy
    assert isinstance(xy, ResultData)
    assert xy.num_fields == 2
    assert len(xy[0].data) == 720
    yz = tensor.yz
    assert isinstance(yz, ResultData)
    assert yz.num_fields == 2
    assert len(yz[0].data) == 720
    xz = tensor.xz
    assert isinstance(xz, ResultData)
    assert xz.num_fields == 2
    assert len(xz[0].data) == 720
    ppal1 = tensor.principal_1
    assert isinstance(ppal1, ResultData)
    assert ppal1.num_fields == 2
    assert len(ppal1[0].data) == 720
    ppal2 = tensor.principal_2
    assert isinstance(ppal2, ResultData)
    assert ppal2.num_fields == 2
    assert len(ppal2[0].data) == 720
    ppal3 = tensor.principal_3
    assert isinstance(ppal3, ResultData)
    assert ppal3.num_fields == 2
    assert len(ppal3[0].data) == 720


class TestCase(unittest.TestCase):
    def test_displacement_elemental_location(self):
        solution = post.load_solution(FILE_PATH)
        self.assertRaises(Exception, solution.displacement, location=locations.elemental)
        try:
            solution.displacement(location=locations.elemental)
        except Exception as e:
            message = "The location must be nodal."
            e2 = Exception(message)
            assert e.args == e2.args
            assert type(e) == type(e2)
            
            
def test_displacement():
    solution = post.load_solution(FILE_PATH)
    vector = solution.displacement()
    print(vector)
    assert vector._operator_name == "U"
    value = vector.vector
    assert value[0].data[3][0] == 9.805953798104982e-06
    assert value[0].data[3][1] == 6.087601335564215e-06
    assert value[0].data[3][2] == -7.841781810225852e-07
    x = vector.x
    assert x[0].data[41] == 9.40664778367545e-07
    y = vector.y
    assert y[0].data[305] == 6.040103203394296e-06
    z = vector.z
    assert z[0].data[548] == -8.479408678360313e-07
    nrm = vector.norm
    assert nrm[0].data[34] ==1.2717854105570665e-06
    

def test_stress():
    solution = post.load_solution(FILE_PATH)
    tensor = solution.stress()
    print(tensor)
    assert tensor._operator_name == "S"
    value = tensor.tensor
    assert value[0].data[3][0] == 12435162.994788285
    assert value[0].data[3][1] == 2.157751263890142e-24
    assert value[0].data[3][2] == -416698.1233175965
    assert value[0].data[3][3] == 5.106217808168383e-09
    assert value[0].data[3][4] == 2.731038343961524e-10
    assert value[0].data[3][5] == 751179.8340354542
    xx = tensor.xx
    assert xx[0].data[41] == 1606559.9330914663
    yy = tensor.yy
    assert yy[0].data[41] == 5.127617333364889e-11
    zz = tensor.zz
    assert zz[0].data[41] == -2856539.4999367334
    xy = tensor.xy
    assert xy[0].data[41] == -3074771.135426731
    yz = tensor.yz
    assert yz[0].data[41] == -75128.58397275927
    xz = tensor.xz
    assert xz[0].data[41] == -1318717.457355396
    ppal1 = tensor.principal_1
    assert ppal1[0].data[41] == 4126211.1436741776
    ppal2 = tensor.principal_2
    assert ppal2[0].data[41] == -1776701.8626811998
    ppal3 = tensor.principal_3
    assert ppal3[0].data[41] == -3599488.8478382444
    vm = tensor.von_mises
    assert vm[0].data[41] == 6994761.422404355



def test_plastic_strain():
    solution = post.load_solution(FILE_PATH)
    tensor = solution.plastic_strain()
    print(tensor)
    assert tensor._operator_name == "EPPL"


def test_elastic_strain():
    solution = post.load_solution(FILE_PATH)
    tensor = solution.elastic_strain()
    print(tensor)
    assert tensor._operator_name == "EPEL"
    value = tensor.tensor
    assert value[0].data[3][0] == 0.00018234479152259739
    assert value[0].data[3][1] == -5.752129670174516e-05
    assert value[0].data[3][2] == -6.555912852208456e-05
    assert value[0].data[3][3] == 1.408180713961465e-19
    assert value[0].data[3][4] == 5.856422410131196e-20
    assert value[0].data[3][5] == 2.8979526592230692e-05
    xx = tensor.xx
    assert xx[0].data[41] == 3.697197923235083e-05
    yy = tensor.yy
    assert yy[0].data[41] ==5.982498289469741e-06
    zz = tensor.zz
    assert zz[0].data[41] == -4.91182636382439e-05
    xy = tensor.xy
    assert xy[0].data[41] == 1.1077761300601659e-20
    yz = tensor.yz
    assert yz[0].data[41] == 1.9300360020838987e-20
    xz = tensor.xz
    assert xz[0].data[41] == -5.087437906548738e-05
    ppal1 = tensor.principal_1
    assert ppal1[0].data[41] == 6.056832330223573e-05
    ppal2 = tensor.principal_2
    assert ppal2[0].data[41] == 5.982498289469729e-06
    ppal3 = tensor.principal_3
    assert ppal3[0].data[41] == -7.271460770812878e-05


def test_temperature():
    solution = post.load_solution(FILE_PATH)
    temp = solution.temperature()
    print(temp)
    assert temp._operator_name == "BFE"
    value = temp.scalar
    assert value.num_fields == 2
    assert len(value[0]) == 720
    assert value[0].data[0] == 22.0
    assert value[0].location == post.locations.nodal
    temp2 = solution.temperature(element_scoping = 2, location = post.locations.elemental)
    value2 = temp2.scalar
    assert value2.num_fields == 1
    assert len(value2[0].data) == 1
    assert value2[0].data[0] == 22.0
    assert value2[0].location == post.locations.elemental