from ansys.dpf.core import locations
import numpy as np
import pytest

from ansys.dpf import post
from ansys.dpf.post import errors
from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.scalar import ComplexScalar, Scalar
from ansys.dpf.post.tensor import ComplexTensor, Tensor
from ansys.dpf.post.vector import ComplexVector, Vector


def test_scalar(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    scalar = solution.structural_temperature()
    assert isinstance(scalar, Scalar)
    assert str(scalar) == (
        "Scalar object.\n\n"
        "Object properties:\n - location   : Nodal\n\n"
        "Structural temperature object."
    )
    value = scalar.scalar
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert value[0].data[0] == 22.0


def test_vector(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    vector = solution.displacement()
    assert isinstance(vector, Vector)
    assert str(vector) == (
        "Vector object.\n\n"
        "Object properties:\n - location   : Nodal\n\n"
        "Displacement object."
    )
    value = vector.vector
    assert isinstance(value, ResultData)
    assert value.num_fields == 1
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 3
    x = vector.x
    assert isinstance(x, ResultData)
    assert x.num_fields == 1
    assert len(x[0].data) > 0
    y = vector.y
    assert isinstance(y, ResultData)
    assert y.num_fields == 1
    assert len(y[0].data) > 0
    z = vector.z
    assert isinstance(z, ResultData)
    assert z.num_fields == 1
    assert len(z[0].data) > 0


def test_tensor(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    tensor = solution.stress()
    assert isinstance(tensor, Tensor)
    assert (
        str(tensor)
        == "Stress Tensor object.\n\nObject properties:\n - location   : Nodal\n"
    )
    value = tensor.tensor
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 6
    xx = tensor.xx
    assert isinstance(xx, ResultData)
    assert xx.num_fields == 2
    assert len(xx[0].data) > 0
    yy = tensor.yy
    assert isinstance(yy, ResultData)
    assert yy.num_fields == 2
    assert len(yy[0].data) > 0
    zz = tensor.zz
    assert isinstance(zz, ResultData)
    assert zz.num_fields == 2
    assert len(zz[0].data) > 0
    xy = tensor.xy
    assert isinstance(xy, ResultData)
    assert xy.num_fields == 2
    assert len(xy[0].data) > 0
    yz = tensor.yz
    assert isinstance(yz, ResultData)
    assert yz.num_fields == 2
    assert len(yz[0].data) > 0
    xz = tensor.xz
    assert isinstance(xz, ResultData)
    assert xz.num_fields == 2
    assert len(xz[0].data) > 0
    ppal1 = tensor.principal_1
    assert isinstance(ppal1, ResultData)
    assert ppal1.num_fields == 2
    assert len(ppal1[0].data) > 0
    ppal2 = tensor.principal_2
    assert isinstance(ppal2, ResultData)
    assert ppal2.num_fields == 2
    assert len(ppal2[0].data) > 0
    ppal3 = tensor.principal_3
    assert isinstance(ppal3, ResultData)
    assert ppal3.num_fields == 2
    assert len(ppal3[0].data) > 0


def test_scalar_complex(complex_model):
    solution = post.load_solution(complex_model)
    scalar = solution.structural_temperature()
    assert isinstance(scalar, ComplexScalar)
    assert str(scalar) == (
        "Complex scalar object.\n\n"
        "Scalar object.\n\n"
        "Object properties:\n - location   : Nodal\n\n"
        "Complex temperature object."
    )
    value = scalar.scalar
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert value[0].data[0] == 22.0
    ampl = scalar.scalar_amplitude
    assert ampl.num_fields == 1
    assert np.isclose(ampl[0].data[0], 31.11269837220809)
    ph = scalar.scalar_at_phase(32.0)
    assert ph.num_fields == 1
    assert np.isclose(ph[0].data[0], 6.9988343023108595)


def test_vector_complex(complex_model):
    solution = post.load_solution(complex_model)
    vector = solution.displacement()
    assert isinstance(vector, ComplexVector)
    assert str(vector) == (
        "Complex vector object.\n\n"
        "Vector object.\n\n"
        "Object properties:\n - location   : Nodal\n"
    )
    value = vector.vector
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 3
    x = vector.x
    assert isinstance(x, ResultData)
    assert x.num_fields == 2
    assert len(x[0].data) > 0
    y = vector.y
    assert isinstance(y, ResultData)
    assert y.num_fields == 2
    assert len(y[0].data) > 0
    z = vector.z
    assert isinstance(z, ResultData)
    assert z.num_fields == 2
    assert len(z[0].data) > 0
    norm = vector.norm
    assert isinstance(z, ResultData)
    assert norm.num_fields == 2
    assert len(norm[0].data) > 0
    value = vector.vector_amplitude
    assert isinstance(value, ResultData)
    assert value.num_fields == 1
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 3
    x = vector.x_amplitude
    assert isinstance(x, ResultData)
    assert x.num_fields == 1
    assert len(x[0].data) > 0
    y = vector.y_amplitude
    assert isinstance(y, ResultData)
    assert y.num_fields == 1
    assert len(y[0].data) > 0
    z = vector.z_amplitude
    assert isinstance(z, ResultData)
    assert z.num_fields == 1
    assert len(z[0].data) > 0
    norm = vector.norm_amplitude
    assert isinstance(z, ResultData)
    assert norm.num_fields == 1
    assert len(norm[0].data) > 0
    ph = vector.vector_at_phase(59.0)
    assert isinstance(ph, ResultData)
    assert ph.num_fields == 1
    assert len(ph[0].data) > 0
    assert len(ph[0].data[3]) == 3
    x_ph = vector.x_at_phase(59.0)
    assert isinstance(x_ph, ResultData)
    assert x_ph.num_fields == 1
    assert len(x_ph[0].data) > 0
    y_ph = vector.y_at_phase(59.0)
    assert isinstance(y_ph, ResultData)
    assert y_ph.num_fields == 1
    assert len(y_ph[0].data) > 0
    z_ph = vector.z_at_phase(59.0)
    assert isinstance(z_ph, ResultData)
    assert z_ph.num_fields == 1
    assert len(z_ph[0].data) > 0
    nrm_ph = vector.norm_at_phase(59.0)
    assert isinstance(nrm_ph, ResultData)
    assert nrm_ph.num_fields == 1
    assert len(nrm_ph[0].data) > 0


def test_tensor_complex(complex_model):
    solution = post.load_solution(complex_model)
    tensor = solution.stress()
    assert isinstance(tensor, ComplexTensor)
    assert str(tensor) == (
        "Complex stress.\n"
        "Complex tensor object.\n\n"
        "Stress Tensor object.\n\n"
        "Object properties:\n - location   : Nodal\n"
    )
    value = tensor.tensor
    assert isinstance(value, ResultData)
    assert value.num_fields == 2
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 6
    xx = tensor.xx
    assert isinstance(xx, ResultData)
    assert xx.num_fields == 2
    assert len(xx[0].data) > 0
    yy = tensor.yy
    assert isinstance(yy, ResultData)
    assert yy.num_fields == 2
    assert len(yy[0].data) > 0
    zz = tensor.zz
    assert isinstance(zz, ResultData)
    assert zz.num_fields == 2
    assert len(zz[0].data) > 0
    xy = tensor.xy
    assert isinstance(xy, ResultData)
    assert xy.num_fields == 2
    assert len(xy[0].data) > 0
    yz = tensor.yz
    assert isinstance(yz, ResultData)
    assert yz.num_fields == 2
    assert len(yz[0].data) > 0
    xz = tensor.xz
    assert isinstance(xz, ResultData)
    assert xz.num_fields == 2
    assert len(xz[0].data) > 0
    ppal1 = tensor.principal_1
    assert isinstance(ppal1, ResultData)
    assert ppal1.num_fields == 2
    assert len(ppal1[0].data) > 0
    ppal2 = tensor.principal_2
    assert isinstance(ppal2, ResultData)
    assert ppal2.num_fields == 2
    assert len(ppal2[0].data) > 0
    ppal3 = tensor.principal_3
    assert isinstance(ppal3, ResultData)
    assert ppal3.num_fields == 2
    assert len(ppal3[0].data) > 0
    value = tensor.tensor_amplitude
    assert isinstance(value, ResultData)
    assert value.num_fields == 1
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 6
    xx = tensor.xx_amplitude
    assert isinstance(xx, ResultData)
    assert xx.num_fields == 1
    assert len(xx[0].data) > 0
    yy = tensor.yy_amplitude
    assert isinstance(yy, ResultData)
    assert yy.num_fields == 1
    assert len(yy[0].data) > 0
    zz = tensor.zz_amplitude
    assert isinstance(zz, ResultData)
    assert zz.num_fields == 1
    assert len(zz[0].data) > 0
    xy = tensor.xy_amplitude
    assert isinstance(xy, ResultData)
    assert xy.num_fields == 1
    assert len(xy[0].data) > 0
    yz = tensor.yz_amplitude
    assert isinstance(yz, ResultData)
    assert yz.num_fields == 1
    assert len(yz[0].data) > 0
    xz = tensor.xz_amplitude
    assert isinstance(xz, ResultData)
    assert xz.num_fields == 1
    assert len(xz[0].data) > 0
    ppal1 = tensor.principal_1_amplitude
    assert isinstance(ppal1, ResultData)
    assert ppal1.num_fields == 1
    assert len(ppal1[0].data) > 0
    ppal2 = tensor.principal_2_amplitude
    assert isinstance(ppal2, ResultData)
    assert ppal2.num_fields == 1
    assert len(ppal2[0].data) > 0
    ppal3 = tensor.principal_3_amplitude
    assert isinstance(ppal3, ResultData)
    assert ppal3.num_fields == 1
    assert len(ppal3[0].data) > 0
    value = tensor.tensor_at_phase(61.0)
    assert isinstance(value, ResultData)
    assert value.num_fields == 1
    assert len(value[0].data) > 0
    assert len(value[0].data[3]) == 6
    xx = tensor.xx_at_phase(61.0)
    assert isinstance(xx, ResultData)
    assert xx.num_fields == 1
    assert len(xx[0].data) > 0
    yy = tensor.yy_at_phase(61.0)
    assert isinstance(yy, ResultData)
    assert yy.num_fields == 1
    assert len(yy[0].data) > 0
    zz = tensor.zz_at_phase(61.0)
    assert isinstance(zz, ResultData)
    assert zz.num_fields == 1
    assert len(zz[0].data) > 0
    xy = tensor.xy_at_phase(61.0)
    assert isinstance(xy, ResultData)
    assert xy.num_fields == 1
    assert len(xy[0].data) > 0
    yz = tensor.yz_at_phase(61.0)
    assert isinstance(yz, ResultData)
    assert yz.num_fields == 1
    assert len(yz[0].data) > 0
    xz = tensor.xz_at_phase(61.0)
    assert isinstance(xz, ResultData)
    assert xz.num_fields == 1
    assert len(xz[0].data) > 0
    ppal1 = tensor.principal_1_at_phase(61.0)
    assert isinstance(ppal1, ResultData)
    assert ppal1.num_fields == 1
    assert len(ppal1[0].data) > 0
    ppal2 = tensor.principal_2_at_phase(61.0)
    assert isinstance(ppal2, ResultData)
    assert ppal2.num_fields == 1
    assert len(ppal2[0].data) > 0
    ppal3 = tensor.principal_3_at_phase(61.0)
    assert isinstance(ppal3, ResultData)
    assert ppal3.num_fields == 1
    assert len(ppal3[0].data) > 0


def test_raise_displacement_elemental_location(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    with pytest.raises(errors.NodalLocationError):
        solution.displacement(location=locations.elemental)


def test_displacement(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    vector = solution.displacement()
    print(vector)
    assert vector._operator_name == "U"
    value = vector.vector
    assert np.isclose(value[0].data[3][0], 9.805953798104982e-06)
    assert np.isclose(value[0].data[3][1], 6.087601335564215e-06)
    assert np.isclose(value[0].data[3][2], -7.841781810225852e-07)
    x = vector.x
    assert np.isclose(x[0].data[41], 9.40664778367545e-07)
    y = vector.y
    assert np.isclose(y[0].data[305], 6.040103203394296e-06)
    z = vector.z
    assert np.isclose(z[0].data[548], -8.479408678360313e-07)
    nrm = vector.norm
    assert np.isclose(nrm[0].data[34], 1.2717854105570665e-06)

    # with dpf.core operator
    from ansys.dpf import core

    op = core.Operator("U")
    # op.inputs.requested_location.connect(core.locations.nodal)
    op.inputs.data_sources.connect(core.DataSources(allkindofcomplexity))
    fc = op.outputs.fields_container()
    assert len(fc) == value.num_fields
    assert fc[0].location == value[0].location
    assert len(fc[0].data) == len(value[0].data)
    assert np.allclose(value[0].data.tolist(), fc[0].data.tolist())
    comp = core.operators.logic.identical_fc()
    comp.inputs.fields_containerA.connect(fc)
    comp.inputs.fields_containerB.connect(value.result_fields_container)
    out = comp.outputs.boolean()
    assert out == True


def test_displacement_complex(complex_model):
    solution = post.load_solution(complex_model)
    vector = solution.displacement()
    print(vector)
    assert vector._operator_name == "U"
    value = vector.vector
    assert np.allclose(
        value[0].data[3].tolist(),
        [2.534371453684853e-09, -5.736467209711275e-10, 6.357980303122968e-11],
    )
    x = vector.x
    assert np.isclose(x[0].data[41], 2.685234654323797e-09)
    y = vector.y
    assert np.isclose(y[0].data[305], -2.442080637069453e-09)
    z = vector.z
    assert np.isclose(z[0].data[548], 1.0919526725085138e-10)
    nrm = vector.norm
    assert np.isclose(nrm[0].data[34], 2.967925671058435e-09)
    value = vector.vector_amplitude
    assert np.allclose(
        value[0].data[3].tolist(),
        [2.5343714759693703e-09, 5.736467469384241e-10, 6.358000469996922e-11],
    )
    x = vector.x_amplitude
    assert np.isclose(x[0].data[41], 2.6852347082946467e-09)
    y = vector.y_amplitude
    assert np.isclose(y[0].data[305], 2.4420806888088805e-09)
    z = vector.z_amplitude
    assert np.isclose(z[0].data[548], 1.0919526860580484e-10)
    nrm = vector.norm_amplitude
    assert np.isclose(nrm[0].data[34], 2.967925756112993e-09)
    value = vector.vector_at_phase(61.0)
    assert np.allclose(
        value[0].data[3].tolist(),
        [1.2283937136871685e-09, -2.7795848616806165e-10, 3.0964159956496574e-11],
    )
    x = vector.x_at_phase(61.0)
    assert np.isclose(x[0].data[41], 1.3013567187124258e-09)
    y = vector.y_at_phase(61.0)
    assert np.isclose(y[0].data[305], -1.183504518054655e-09)
    z = vector.z_at_phase(61.0)
    assert np.isclose(z[0].data[548], 5.292387083515219e-11)
    nrm = vector.norm_at_phase(61.0)
    assert np.isclose(nrm[0].data[34], 1.438258083761136e-09)


def test_stress(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    tensor = solution.stress()
    assert tensor._operator_name == "S"
    value = tensor.tensor
    assert value[0].data[3][0]
    assert value[0].data[3][1]
    assert value[0].data[3][2]
    assert value[0].data[3][3]
    assert value[0].data[3][4]
    assert value[0].data[3][5]
    xx = tensor.xx
    assert xx[0].data[41]
    yy = tensor.yy
    assert yy[0].data[41]
    zz = tensor.zz
    assert zz[0].data[41]
    xy = tensor.xy
    assert xy[0].data[41]
    yz = tensor.yz
    assert yz[0].data[41]
    xz = tensor.xz
    assert xz[0].data[41]
    ppal1 = tensor.principal_1
    assert ppal1[0].data[41]
    ppal2 = tensor.principal_2
    assert ppal2[0].data[41]
    ppal3 = tensor.principal_3
    assert ppal3[0].data[41]
    vm = tensor.von_mises
    assert vm[0].data[41]

    # with dpf.core operator
    from ansys.dpf import core

    op = core.Operator("S")
    op.inputs.requested_location.connect(core.locations.nodal)
    op.inputs.data_sources.connect(core.DataSources(allkindofcomplexity))
    fc = op.outputs.fields_container()
    assert len(fc) == value.num_fields
    assert fc[0].location == value[0].location
    assert len(fc[0].data) == len(value[0].data)
    assert np.allclose(value[0].data.tolist(), fc[0].data.tolist())
    comp = core.operators.logic.identical_fc()
    comp.inputs.fields_containerA.connect(fc)
    comp.inputs.fields_containerB.connect(value.result_fields_container)
    out = comp.outputs.boolean()
    assert out is True


def test_stress_complex(complex_model):
    solution = post.load_solution(complex_model)
    tensor = solution.stress()
    print(tensor)
    assert tensor._operator_name == "S"
    value = tensor.tensor
    assert np.allclose(
        value[0].data[3].tolist(),
        [
            -1894.3998413085938,
            -99533.1953125,
            -216.0846405029297,
            -15840.79736328125,
            548.1216735839844,
            -3538.7244873046875,
        ],
    )
    xx = tensor.xx
    yy = tensor.yy
    zz = tensor.zz
    xy = tensor.xy
    yz = tensor.yz
    xz = tensor.xz
    ppal1 = tensor.principal_1
    ppal2 = tensor.principal_2
    ppal3 = tensor.principal_3
    vm = tensor.von_mises
    assert vm.num_fields == 2
    assert len(vm[0].data) == 4802
    assert np.isclose(xx[0].data[41], -41958.8046875)
    assert np.isclose(yy[0].data[41], -75517.203125)
    assert np.isclose(zz[0].data[41], -793.2758483886719)
    assert np.isclose(xy[0].data[41], -50654.06640625)
    assert np.isclose(yz[0].data[41], -3254.9710693359375)
    assert np.isclose(xz[0].data[41], 4329.180419921875)
    assert np.isclose(ppal1[0].data[41], 2795.2100794761645)
    assert np.isclose(ppal2[0].data[41], -8965.582783669306)
    assert np.isclose(ppal3[0].data[41], -112098.91095669553)
    assert np.isclose(vm[0].data[41], 109488.4895262726)
    xx = tensor.xx_amplitude
    yy = tensor.yy_amplitude
    zz = tensor.zz_amplitude
    xy = tensor.xy_amplitude
    yz = tensor.yz_amplitude
    xz = tensor.xz_amplitude
    ppal1 = tensor.principal_1_amplitude
    ppal2 = tensor.principal_2_amplitude
    ppal3 = tensor.principal_3_amplitude
    vm = tensor.von_mises_amplitude
    assert vm.num_fields == 1
    assert len(vm[0].data) > 0
    assert np.isclose(xx[0].data[41], 41958.80541595527)
    assert np.isclose(yy[0].data[41], 75517.21861853577)
    assert np.isclose(zz[0].data[41], 793.3104230072076)
    assert np.isclose(xy[0].data[41], 50654.111124516836)
    assert np.isclose(yz[0].data[41], 3255.0650347243854)
    assert np.isclose(xz[0].data[41], 4329.302877644564)
    assert np.isclose(ppal1[0].data[41], 2796.0674775873003)
    assert np.isclose(ppal2[0].data[41], 8965.60427429267)
    assert np.isclose(ppal3[0].data[41], 112098.95413919646)
    assert np.isclose(vm[0].data[41], 109488.58588814907)
    xx = tensor.xx_at_phase(9.0)
    yy = tensor.yy_at_phase(9.0)
    zz = tensor.zz_at_phase(9.0)
    xy = tensor.xy_at_phase(9.0)
    yz = tensor.yz_at_phase(9.0)
    xz = tensor.xz_at_phase(9.0)
    ppal1 = tensor.principal_1_at_phase(9.0)
    ppal2 = tensor.principal_2_at_phase(9.0)
    ppal3 = tensor.principal_3_at_phase(9.0)
    vm = tensor.von_mises_at_phase(9.0)
    assert vm.num_fields == 1
    assert len(vm[0].data) > 0
    assert np.isclose(xx[0].data[41], -41440.999079450754)
    assert np.isclose(yy[0].data[41], -74579.8936585364)
    assert np.isclose(zz[0].data[41], -784.6679315714965)
    assert np.isclose(xy[0].data[41], -50019.90154956536)
    assert np.isclose(yz[0].data[41], -3218.7660576056496)
    assert np.isclose(xz[0].data[41], 4280.974878496336)
    assert np.isclose(ppal1[0].data[41], 2749.9651388433817)
    assert np.isclose(ppal2[0].data[41], -8852.130711811371)
    assert np.isclose(ppal3[0].data[41], -110703.39509659067)
    assert np.isclose(vm[0].data[41], 108117.78055484823)


def test_plastic_strain(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    tensor = solution.plastic_strain()
    print(tensor)
    assert tensor._operator_name == "EPPL"


def test_plastic_strain_complex(complex_model):
    solution = post.load_solution(complex_model)
    tensor = solution.plastic_strain()
    print(tensor)
    assert tensor._operator_name == "EPPL"


def test_elastic_strain(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    tensor = solution.elastic_strain()
    print(tensor)
    assert tensor._operator_name == "EPEL"
    value = tensor.tensor
    assert value[0].data[3][0]
    assert value[0].data[3][1]
    assert value[0].data[3][2]
    assert value[0].data[3][3]
    assert value[0].data[3][4]
    assert value[0].data[3][5]
    xx = tensor.xx
    assert xx[0].data[41]
    yy = tensor.yy
    assert yy[0].data[41]
    zz = tensor.zz
    assert zz[0].data[41]
    xy = tensor.xy
    assert xy[0].data[41]
    yz = tensor.yz
    assert yz[0].data[41]
    xz = tensor.xz
    assert xz[0].data[41]
    ppal1 = tensor.principal_1
    assert ppal1[0].data[41]
    ppal2 = tensor.principal_2
    assert ppal2[0].data[41]
    ppal3 = tensor.principal_3
    assert ppal3[0].data[41]

    # with dpf.core operator
    from ansys.dpf import core

    op = core.Operator("EPEL")
    op.inputs.requested_location.connect(core.locations.nodal)
    op.inputs.data_sources.connect(core.DataSources(allkindofcomplexity))
    fc = op.outputs.fields_container()
    assert len(fc) == value.num_fields
    assert fc[0].location == value[0].location
    assert len(fc[0].data) == len(value[0].data)
    assert np.allclose(value[0].data.tolist(), fc[0].data.tolist())
    comp = core.operators.logic.identical_fc()
    comp.inputs.fields_containerA.connect(fc)
    comp.inputs.fields_containerB.connect(value.result_fields_container)
    out = comp.outputs.boolean()
    assert out is True


def test_elastic_strain_complex(complex_model):
    solution = post.load_solution(complex_model)
    tensor = solution.elastic_strain()
    print(tensor)
    assert tensor._operator_name == "EPEL"
    value = tensor.tensor
    assert np.allclose(
        value[0].data[3].tolist(),
        [
            3.031909585615722e-07,
            -7.12252500534305e-07,
            3.211454924212376e-07,
            -1.400326468115054e-07,
            2.5393885882962763e-09,
            -2.6922432283527087e-08,
        ],
    )
    assert tensor.xx
    assert tensor.yy
    assert tensor.zz
    assert tensor.xy
    assert tensor.yz
    assert tensor.xz
    assert tensor.principal_1
    assert tensor.principal_2
    assert tensor.principal_3
    assert tensor.xx_amplitude
    assert tensor.yy_amplitude
    assert tensor.zz_amplitude
    assert tensor.xy_amplitude
    assert tensor.yz_amplitude
    assert tensor.xz_amplitude
    assert tensor.principal_1_amplitude
    assert tensor.principal_2_amplitude
    assert tensor.principal_3_amplitude
    assert tensor.xx_at_phase(14.0)
    assert tensor.yy_at_phase(14.0)
    assert tensor.zz_at_phase(14.0)
    assert tensor.xy_at_phase(14.0)
    assert tensor.yz_at_phase(14.0)
    assert tensor.xz_at_phase(14.0)
    assert tensor.principal_1_at_phase(14.0)
    assert tensor.principal_2_at_phase(14.0)
    assert tensor.principal_3_at_phase(14.0)


def test_temperature(allkindofcomplexity):
    solution = post.load_solution(allkindofcomplexity)
    temp = solution.structural_temperature()
    print(temp)
    assert temp._operator_name == "BFE"
    value = temp.scalar
    assert value.num_fields == 2
    assert value[0].data[0] == 22.0
    assert value[0].location == post.locations.nodal
    temp2 = solution.structural_temperature(
        element_scoping=2, location=post.locations.elemental
    )
    value2 = temp2.scalar
    assert value2.num_fields == 1
    assert len(value2[0].data) == 1
    assert value2[0].data[0] == 22.0
    assert value2[0].location == post.locations.elemental

    # with dpf.core operator
    from ansys.dpf import core

    op = core.Operator("BFE")
    op.inputs.requested_location.connect(core.locations.nodal)
    op.inputs.data_sources.connect(core.DataSources(allkindofcomplexity))
    fc = op.outputs.fields_container()
    assert len(fc) == value.num_fields
    assert fc[0].location == value[0].location
    assert len(fc[0].data) == len(value[0].data)
    assert np.allclose(value[0].data.tolist(), fc[0].data.tolist())
    comp = core.operators.logic.identical_fc()
    comp.inputs.fields_containerA.connect(fc)
    comp.inputs.fields_containerB.connect(value.result_fields_container)
    out = comp.outputs.boolean()
    assert out is True


def test_temperature_complex(complex_model):
    solution = post.load_solution(complex_model)
    temp = solution.structural_temperature()
    print(temp)
    assert temp._operator_name == "BFE"
    value = temp.scalar
    assert value.num_fields == 2
    assert len(value[0]) > 0
    assert value[0].data[0] == 22.0
    assert value[0].location == post.locations.nodal
    value = temp.scalar_amplitude
    assert value.num_fields == 1
    assert len(value[0]) > 0
    assert np.isclose(value[0].data[0], 31.11269837220809)
    assert value[0].location == post.locations.nodal
    value = temp.scalar_at_phase(21.0)
    assert value.num_fields == 1
    assert len(value[0]) > 0
    assert np.isclose(value[0].data[0], 12.654674492941831)
    assert value[0].location == post.locations.nodal
    temp2 = solution.structural_temperature(
        element_scoping=2, location=post.locations.elemental
    )
    value2 = temp2.scalar
    assert value2.num_fields == 2
    assert len(value2[0].data) == 1
    assert value2[0].data[0] == 22.0
    assert value2[0].location == post.locations.elemental
    value2 = temp2.scalar_amplitude
    assert value2.num_fields == 1
    assert len(value2[0].data) == 1
    assert np.isclose(value2[0].data[0], 31.11269837220809)
    assert value2[0].location == post.locations.elemental
    value2 = temp2.scalar_at_phase(21.0)
    assert value2.num_fields == 1
    assert len(value2[0].data) == 1
    assert np.isclose(value2[0].data[0], 12.654674492941835)
    assert value2[0].location == post.locations.elemental
