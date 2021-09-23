"""
.. _ref_harmonic_analysis:

ANSYS DPF-Post: Harmonic Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows how post-process a harmonic analysis result file
using API of the DPF-POST module.
"""

###############################################################################
# Get started
# ~~~~~~~~~~~
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get the solution object
# ~~~~~~~~~~~~~~~~~~~~~~~
# The following file is the result of a harmonic analysis computed
# using Ansys Mechanical.
#
# Here we load the solution
solution = post.load_solution(examples.complex_rst)
print(solution)

###############################################################################
# Get result objects
# ~~~~~~~~~~~~~~~~~~

###############################################################################
# Get a displacement result
# =========================

###############################################################################
# **Get the result**: it will contain a field for real values and a
# field for imaginary values.
disp_result = solution.displacement()
disp = disp_result.vector

###############################################################################
# **Check the number of fields**
disp.num_fields

###############################################################################
# **Get data from a field**
disp.get_data_at_field(0)

###############################################################################
# **Get maximum data value over all fields**
disp.max_data

###############################################################################
# **Get minimum data value over all fields**
disp.min_data

###############################################################################
# **Get maximum data value over a targeted field**
disp.get_max_data_at_field(0)

###############################################################################
# **Get minimum data value over all fields**
disp.get_min_data_at_field(0)

###############################################################################
# Get a stress result and deals with amplitude
# ============================================

###############################################################################
# **Get the result**
stress_result = solution.stress()

###############################################################################
# **Check if the support has complex frequencies**
stress_result.has_complex_frequencies()

###############################################################################
# **Get the tensor result**
stress = stress_result.tensor
# shell and solid elements are in distinct fields.
# so we have four fields : the solid-real one, the solid-imaginary one,
# the shell-real one and the shell_imaginary one.
stress.num_fields

###############################################################################
# **Get the shell field**
shell_field = stress[0]
shell_field.shell_layers

###############################################################################
# **Get the solid field field**
solid_field = stress[1]

###############################################################################
# **Plot the amplitude contour**
amplitude = stress_result.tensor_amplitude
amplitude.plot_contour()

###############################################################################
# Get an elastic_strain result and deals with phase
# =================================================

###############################################################################
# **Get the result**: it will contain a field for real values and a
# field for imaginary values.
elastic_strain_result = solution.elastic_strain()
elastic_strain = elastic_strain_result.tensor
# shell and solid elements are in distinct fields.
elastic_strain.num_fields

###############################################################################
# **Deal with phase**: phase unit is degree, phase must be a float.
es_at_phase = elastic_strain_result.tensor_at_phase(39.0)
es_at_phase.max_data
es_at_phase.num_fields
real_field = elastic_strain_result.tensor_at_phase(0.0)
img_field = elastic_strain_result.tensor_at_phase(90.0)

###############################################################################
# **It is also possible to deal with plastic_strain and temperature this way.**
# The result file must contain those results.
