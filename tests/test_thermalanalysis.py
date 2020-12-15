from ansys.dpf.post.temperature import Temperature
from ansys.dpf.post.common import _PhysicsType
from ansys.dpf import post
import numpy as np
import unittest
import pytest

filepath_steady = "d:/rst/rth_steady.rth"
filepath_transient = "d:/rst/rth_transient.rth"


def test_thermal_steadystate():
    solution = post.load_solution(filepath_steady)
    assert solution._model.metadata.result_info.physics_type == _PhysicsType.thermal
    temp = solution.temperature()
    assert isinstance(temp, Temperature)
    s = temp.scalar
    assert s.num_fields == 1
    assert np.isclose(s[0].data[23], 29.6247641917003)

to_return = "node scoping and element scoping returns the same"
def test_steadystate_nodscoping():
    solution = post.load_solution(filepath_steady)
    temp = solution.temperature(node_scoping = [2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    #assert np.allclose(s[0].data.tolist(), [41.8243869 , 40.29943406])
    raise Exception(to_return)


def test_steadystate_elemscoping():
    solution = post.load_solution(filepath_steady)
    temp = solution.temperature(element_scoping = [2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    #assert np.allclose(s[0].data.tolist(), [41.8243869 , 40.29943406])
    raise Exception(to_return)


def test_steadystate_nodlocation():
    solution = post.load_solution(filepath_steady)
    temp = solution.temperature(location = post.locations.nodal)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal

class TestCase1(unittest.TestCase):
    
    @pytest.fixture(autouse=True)
    def set_filepath(self, allkindofcomplexity):
        self._filepath = filepath_steady
        
    def test_steadystate_elemlocation(self):
        solution = post.load_solution(self._filepath)
        temp = solution.temperature(location = post.locations.elemental)
        self.assertRaises(Exception, s = temp.scalar)
    
    def test_steadystate_elemnodallocation(self):
        solution = post.load_solution(self._filepath)
        temp = solution.temperature(location = post.locations.elemental_nodal)
        self.assertRaises(Exception, s = temp.scalar)


def test_steadystate_timescoping():
    solution = post.load_solution(filepath_steady)
    temp = solution.temperature(time_scoping = 1)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 30.387240610202973)


def test_steadystate_time():
    solution = post.load_solution(filepath_steady)
    temp = solution.temperature(time = 1.)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 30.387240610202973)


def test_steadystate_set():
    solution = post.load_solution(filepath_steady)
    temp = solution.temperature(set = 1)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 30.387240610202973)


def test_thermal_transient():
    solution = post.load_solution(filepath_transient)
    assert solution._model.metadata.result_info.physics_type == _PhysicsType.thermal
    temp = solution.temperature()
    assert isinstance(temp, Temperature)
    s = temp.scalar
    assert s.num_fields == 1
    assert np.isclose(s[0].data[23], 22.159122145252788)


def test_transient_nodscoping():
    solution = post.load_solution(filepath_transient)
    temp = solution.temperature(node_scoping = [2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    #assert np.allclose(s[0].data.tolist(), [27.01872925, 25.61222966])
    raise Exception(to_return)


def test_transient_elemscoping():
    solution = post.load_solution(filepath_transient)
    temp = solution.temperature(node_scoping = [2, 45])
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 2
    #assert np.allclose(s[0].data.tolist(), [27.01872925, 25.61222966])
    raise Exception(to_return)


def test_transient_nodlocation():
    solution = post.load_solution(filepath_transient)
    temp = solution.temperature(location = post.locations.nodal)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal


class TestCase2(unittest.TestCase):
    
    @pytest.fixture(autouse=True)
    def set_filepath(self, allkindofcomplexity):
        self._filepath = filepath_transient
        
    def test_steadystate_elemlocation(self):
        solution = post.load_solution(self._filepath)
        temp = solution.temperature(location = post.locations.elemental)
        self.assertRaises(Exception, s = temp.scalar)
    
    def test_steadystate_elemnodallocation(self):
        solution = post.load_solution(self._filepath)
        temp = solution.temperature(location = post.locations.elemental_nodal)
        self.assertRaises(Exception, s = temp.scalar)


def test_transient_timescoping():
    solution = post.load_solution(filepath_transient)
    temp = solution.temperature(time_scoping = [2, 15])
    s = temp.scalar
    assert s.num_fields == 2
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 21.999999999992536)


def test_transient_time():
    solution = post.load_solution(filepath_transient)
    temp = solution.temperature(set = 2)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 21.999999999992536)


def test_transient_set():
    solution = post.load_solution(filepath_transient)
    temp = solution.temperature(time = 0.0223)
    s = temp.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[24], 21.999999999992323)