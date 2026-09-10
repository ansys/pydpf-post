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
.. _ref_complex_results:

Access complex results from a modal or harmonic analysis
--------------------------------------------------------
This example shows how to use the legacy PyDPF-Post API to access complex results
from a modal or harmonic analysis.
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
# Get the ``Solution`` object. This example loads a file with complex results.

solution = post.load_solution(examples.complex_rst)
solution.has_complex_result()

###############################################################################
# Get displacement result
# -----------------------
# The displacement result contain a field for real values and a field for
# imaginary values.

disp_result = solution.displacement()

###############################################################################
# Check for support of complex frequencies
# ----------------------------------------

disp_result.has_complex_frequencies()

###############################################################################
# Compute result
# --------------
disp = disp_result.vector
disp.num_fields

###############################################################################
# Define phase
# ------------
# The phase value must be a float. The phase unit is degrees.

phase = 39.0
disp_at_phase = disp_result.vector_at_phase(phase)
print(f"Maximum displacement at phase {phase}: ", disp_at_phase.max_data)
print(f"There are {disp_at_phase.num_fields} fields")
real_field = disp_result.vector_at_phase(0.0)
img_field = disp_result.vector_at_phase(90.0)

real_field

###############################################################################
# Get amplitude
# -------------

disp_ampl = disp_result.vector_amplitude
disp_ampl.num_fields
disp_ampl.max_data
