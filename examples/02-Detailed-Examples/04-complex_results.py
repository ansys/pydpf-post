"""
.. _ref_complex_results:

Compute Complex Results from Modal or Harmonic Analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows how to access complex results using the DPF-POST
module API.
"""

###############################################################################
# **Get started**
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# **Get the solution object**
solution = post.load_solution(examples.complex_rst)
solution.has_complex_result()

###############################################################################
# **Get displacement result**
#
# This result will contain a field for real values and a field for
# imaginary values)
disp_result = solution.displacement()

###############################################################################
# **Check if the support has complex frequencies**
disp_result.has_complex_frequencies()

###############################################################################
# **Compute the result**
disp = disp_result.vector
disp.num_fields

###############################################################################
# **Deals with phase** (phase unit is degree, phase must be a float)
phase = 39.0
disp_at_phase = disp_result.vector_at_phase(phase)
print(f"Maximum displacement at phase {phase}:", disp_at_phase.max_data)
print(f"There are {disp_at_phase.num_fields} fields")
real_field = disp_result.vector_at_phase(0.0)
img_field = disp_result.vector_at_phase(90.0)

real_field

###############################################################################
# **Get amplitude**
disp_ampl = disp_result.vector_amplitude
disp_ampl.num_fields
disp_ampl.max_data
