"""
.. _ref_complex_results:

Complex results from a modal or harmonic analysis
-------------------------------------------------
This example shows how you can access complex results from a modal or
harmonic analysis.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Solution`` object
# -----------------------
# Get the ``Solution`` object.

solution = post.load_solution(examples.complex_rst)
solution.has_complex_result()

###############################################################################
# Get displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~
# Get the displacement ``Result`` object. It contain a field for real values
# and a field for imaginary values.

disp_result = solution.displacement()

###############################################################################
# Check if support has complex frequencies
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check if the support has complex frequencies.

disp_result.has_complex_frequencies()

###############################################################################
# **Compute the result**
disp = disp_result.vector
disp.num_fields

###############################################################################
# Define phase
# ~~~~~~~~~~~~
# Define the phase. The phase value must be a float. The phase unit is degrees.

phase = 39.0
disp_at_phase = disp_result.vector_at_phase(phase)
print(f"Maximum displacement at phase {phase}:", disp_at_phase.max_data)
print(f"There are {disp_at_phase.num_fields} fields")
real_field = disp_result.vector_at_phase(0.0)
img_field = disp_result.vector_at_phase(90.0)

real_field

###############################################################################
# Get amplitude
# ~~~~~~~~~~~~~
# Get the amplitude.

disp_ampl = disp_result.vector_amplitude
disp_ampl.num_fields
disp_ampl.max_data
