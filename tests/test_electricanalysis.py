from ansys.dpf.post.electric_results import ElectricField, ElectricPotential
from ansys.dpf.post.common import _PhysicsType
from ansys.dpf import post
import numpy as np
import unittest
import pytest
from ansys.dpf.post import errors as dpf_errors
from ansys.dpf.core import errors as core_errors


def test_electricfield(rth_electric):
    solution = post.load_solution(rth_electric)
    assert solution._model.metadata.result_info.physics_type == _PhysicsType.thermal
    ef = solution.electric_field()
    assert isinstance(ef, ElectricField)
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data[20]) == 3
    assert np.isclose(s[0].data[23][1], 19.562952041625977)


def test_electricfield_nodscoping(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(node_scoping = [2])
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data) == 1
    assert len(s[0].data[0]) == 3
    assert np.allclose(s[0].data[0].tolist(), [5.25223311e-14, 1.95629520e+01, 2.82945325e-14])
    ef = solution.electric_field(location = post.locations.elemental, node_scoping = [2])
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.elemental
    assert len(s[0].data) == 8
    assert len(s[0].data[0]) == 3
    assert np.allclose(s[0].data[0].tolist(), [-3.41948692e-14,  1.95629520e+01,  7.77156117e-15])
    ef = solution.electric_field(location = post.locations.elemental_nodal, node_scoping = [2])
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.elemental_nodal
    assert len(s[0].data) == 8
    assert len(s[0].data[0]) == 3
    assert np.allclose(s[0].data.tolist(), [ 2.63128894e-11,  1.95629520e+01,  2.62733394e-11])


@pytest.mark.skipif(True, reason="element scoping not available with electrical results.")
def test_electricfield_elemscoping(rth_electric):
    raise Exception("Element scoping on electric_field does not work.")
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(element_scoping = [2])
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data) == 20
    assert len(s[0].data[0]) == 3
    #assert np.isclose(s[0].data[0].tolist(), [2.63128894e-11, 1.95629520e+01, 2.62733394e-11])
    ef = solution.electric_field(location = post.locations.elemental, element_scoping = [2])
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.elemental
    assert len(s[0].data) == 3
    #assert np.isclose(s[0].data.tolist(), [-3.41948692e-14,  1.95629520e+01,  7.77156117e-15])
    ef = solution.electric_field(location = post.locations.elemental_nodal, element_scoping = [2])
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.elemental_nodal
    assert len(s[0].data) == 8
    assert len(s[0].data[0]) == 3
    #assert np.isclose(s[0].data.tolist(), [-3.41948692e-14,  1.95629520e+01,  7.77156117e-15])


def test_electricfield_nodlocation(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field()
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal


def test_electricfield_elemlocation(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(location = post.locations.elemental)
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.elemental


def test_electricfield_elemnodlocation(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(location = post.locations.elemental_nodal)
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.elemental_nodal


def test_electricfield_timescoping(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(time_scoping = 1)
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data[20]) == 3
    assert np.isclose(s[0].data[23][1], 19.562952041625977)


def test_electricfield_time(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(time = 1.)
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data[20]) == 3
    assert np.isclose(s[0].data[23][1], 19.562952041625977)


def test_electricfield_set(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_field(set = 1)
    s = ef.vector
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data[20]) == 3
    assert np.isclose(s[0].data[23][1], 19.562952041625977)



def test_electricpotential(rth_electric):
    solution = post.load_solution(rth_electric)
    assert solution._model.metadata.result_info.physics_type == _PhysicsType.thermal
    ef = solution.electric_potential()
    assert isinstance(ef, ElectricPotential)
    s = ef.scalar
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data) == 4125
    assert np.isclose(s[0].data[23], 0.09781476007338061)
    

to_return = "node scoping and element scoping returns the same"
def test_electricpotential_nodscoping(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_potential(node_scoping = [2])
    s = ef.scalar
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data) == 1
    assert np.isclose(s[0].data[0], 0.048907380036668786)


@pytest.mark.skipif(True, reason="element scoping not available with electrical results.")
def test_electricpotential_elemscoping(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_potential(node_scoping = [2])
    s = ef.scalar
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal
    assert len(s[0].data) == 1
    #assert np.isclose(s[0].data[0], 0.02445369)
    raise Exception(to_return)


def test_electricpotential_nodlocation(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_potential(location = post.locations.nodal)
    s = ef.scalar
    assert s.num_fields == 1
    assert s[0].location == post.locations.nodal

  
def test_electricpotential_elemlocation(rth_electric):
    solution = post.load_solution(rth_electric)
    with pytest.raises(dpf_errors.NodalLocationError):
        pot = solution.electric_potential(location = post.locations.elemental)
        

def test_electricpotential_elemnodallocation(rth_electric):
    solution = post.load_solution(rth_electric)
    with pytest.raises(dpf_errors.NodalLocationError):
        pot = solution.electric_potential(location = post.locations.elemental_nodal)


def test_electricpotential_timescoping(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_potential(time_scoping = [1])
    s = ef.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[0], 0.07336107005500624)


def test_electricpotential_time(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_potential(set = 1)
    s = ef.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[0], 0.07336107005500624)


def test_electricpotential_set(rth_electric):
    solution = post.load_solution(rth_electric)
    ef = solution.electric_potential(time = 1.)
    s = ef.scalar
    assert s.num_fields == 1
    assert len(s[0].data) == 4125
    assert s[0].location == post.locations.nodal
    assert np.isclose(s[0].data[0], 0.07336107005500624)
    
    
if __name__ == "__main__":
    rth_electric = "d:/AnsysDev/DPF-Core/ansys/dpf/core/examples/rth/rth_electric.rth"
    test_electricfield(rth_electric)
    