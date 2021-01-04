"""
.. _ref_complex_results:

ANSYS DPF Computes complex result with the POST Api
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows how to access the complex API 
of the POST module.
"""

###############################################################################
# **Get started**
from ansys.dpf import post


###############################################################################
# **Get the solution object**
import os
path = os.getcwd()
path += "/../../tests/testfiles/complex/fileComplex.rst"
solution = post.load_solution(path)
solution.has_complex_result()

###############################################################################
# **Get displacement result** (it will contain a field for real values and a 
# field for imaginary values)
disp_result = solution.displacement()
disp = disp_result.vector
disp.is_complex_result()
disp.num_fields

###############################################################################
# **Deals with phase** (phase unit is degree, phase must be a float)
disp_at_phase = disp_result.vector_at_phase(39.)
disp_at_phase.max_data
disp_at_phase.num_fields
real_field = disp_result.vector_at_phase(0.)
img_field = disp_result.vector_at_phase(90.)

###############################################################################
# **Get amplitude**
disp_ampl = disp_result.vector_amplitude
disp_ampl.max_data
disp_ampl.num_fields

