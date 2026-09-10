"""
.. _ref_modal_analysis:

Modal analysis
==============
This example shows how you can postprocess the result file for a modal analysis
using DPF-Post.
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
# Get the ``Solution`` object. This example loads a result file for a modal analysis
# computed in Ansys Mechanical.

example_path = examples.download_all_kinds_of_complexity_modal()

solution = post.load_solution(example_path)
print(solution)


###############################################################################
# Get ``Result`` objects
# ----------------------

###############################################################################
# Get displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~
# Get the displacement ``Result`` object. It contains a field for real values
# and a field for imaginary values.

disp_result = solution.displacement()
disp = disp_result.vector

###############################################################################
# Check number of fields
# ~~~~~~~~~~~~~~~~~~~~~~
# Check the number of fields.

disp.num_fields

###############################################################################
# **Get data from a field**

disp.get_data_at_field(0)

###############################################################################
# Get maximum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the maximum data value over all fields.

disp.max_data

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the minimum data value over all fields.

disp.min_data

###############################################################################
# Get maximum data value over targeted field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the maximum data value over a targeted field.

disp.get_max_data_at_field(0)

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the minimum data value over all fields.

disp.get_min_data_at_field(0)

###############################################################################
# Get stress result
# -----------------
# Get a stress result that deals with amplitude. It contains a field for real
# values and a field for imaginary values.

stress_result = solution.stress()

###############################################################################
# Check if support has complex frequencies
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check if the support has complex frequencies.

stress_result.has_complex_frequencies()

###############################################################################
# Get tensor result
# ~~~~~~~~~~~~~~~~~
# Get the tensor result.

stress = stress_result.tensor
stress.num_fields

###############################################################################
# Get shell field
# ~~~~~~~~~~~~~~~
# Get the shell field.

shell_field = stress[0]
shell_field.shell_layers

###############################################################################
# Get solid field
# ~~~~~~~~~~~~~~~
# Get the solid field.

solid_field = stress[1]

###############################################################################
# Plot amplitude contour
# ~~~~~~~~~~~~~~~~~~~~~~
# Plot the amplitude contour.

amplitude = stress_result.tensor_amplitude
stress.plot_contour()

###############################################################################
# Get elastic strain result
# =========================
# Get an elastic strain result that deals with phase. It contains a field for
# real values and a field for imaginary values.

elastic_strain_result = solution.elastic_strain()
elastic_strain = elastic_strain_result.tensor

# Check number of fields
# ~~~~~~~~~~~~~~~~~~~~~~
# Check the number of shell and solid elements in distinct fields.

elastic_strain.num_fields

###############################################################################
# If the result file contains results, you can use this method
# to get the elastic strain result.

print(solution.plastic_strain())
