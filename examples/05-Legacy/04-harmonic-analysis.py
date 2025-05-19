# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
.. _ref_harmonic_analysis:

Postprocess the result file for a harmonic analysis
===================================================
This example shows how to use the legacy PyDPF-Post API to postprocess a result
file for a harmonic analysis.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Solution`` object
# -----------------------
# Get the ``Solution`` object. This example loads a result file for an harmonic
# analysis computed in Ansys Mechanical.

example_path = examples.download_all_kinds_of_complexity()

solution = post.load_solution(examples.complex_rst)
print(solution)

###############################################################################
# Get ``Result`` objects
# ----------------------

###############################################################################
# Get displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~
# The displacement result contains a field for real values and a field for
# imaginary values.

disp_result = solution.displacement()
disp = disp_result.vector

###############################################################################
# Get number of fields
# ~~~~~~~~~~~~~~~~~~~~

disp.num_fields

###############################################################################
# Get data from field
# ~~~~~~~~~~~~~~~~~~~
disp.get_data_at_field(0)

###############################################################################
# Get maximum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.max_data

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.min_data

###############################################################################
# Get maximum data value over targeted field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.get_max_data_at_field(0)

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.get_min_data_at_field(0)

###############################################################################
# Get stress result
# ~~~~~~~~~~~~~~~~~
# Get a stress result that deals with amplitude. It contains a field for real
# values and a field for imaginary values.

stress_result = solution.stress()

###############################################################################
# Check if support has complex frequencies
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stress_result.has_complex_frequencies()

###############################################################################
# Get tensor result
# ~~~~~~~~~~~~~~~~~

stress = stress_result.tensor

###############################################################################
# Get number of fields
# ~~~~~~~~~~~~~~~~~~~~~~
# Get the number of shell and solid elements in distinct fields. Because shell and
# solid elements are in distinct fields, you get four fields: the
# solid-real one, the solid-imaginary one, the shell-real one, and the
# shell-imaginary one.

stress.num_fields

###############################################################################
# Get shell field
# ~~~~~~~~~~~~~~~

shell_field = stress[0]
shell_field.shell_layers

###############################################################################
# Get solid field
# ~~~~~~~~~~~~~~~~

solid_field = stress[1]

###############################################################################
# Plot amplitude contour
# ~~~~~~~~~~~~~~~~~~~~~~

amplitude = stress_result.tensor_amplitude
amplitude.plot_contour()

###############################################################################
# Get elastic strain result
# -------------------------
# Get an elastic strain result that deals with phase. It contains a field for
# real values and a field for imaginary values.

elastic_strain_result = solution.elastic_strain()
elastic_strain = elastic_strain_result.tensor
# shell and solid elements are in distinct fields.
elastic_strain.num_fields

###############################################################################
# Define phase
# ~~~~~~~~~~~~
# The phase must be a float value. The unit is degrees.

es_at_phase = elastic_strain_result.tensor_at_phase(39.0)
es_at_phase.max_data
es_at_phase.num_fields
real_field = elastic_strain_result.tensor_at_phase(0.0)
img_field = elastic_strain_result.tensor_at_phase(90.0)

###############################################################################
# If the result file contains results, use this method
# to get the elastic strain result.

print(solution.elastic_strain())
