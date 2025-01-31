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
.. _ref_result_on_path:

Request result on a specific path
=================================
This example shows how you can request a result on a
specific path of coordinates.

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

solution = post.load_solution(examples.static_rst)

###############################################################################
# Create coordinates array
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Create the coordinates array to request the result on.

coordinates = [[0.024, 0.03, 0.003]]
for i in range(1, 51):
    coord_copy = coordinates[0].copy()
    coord_copy[1] = coord_copy[0] + i * 0.001
    coordinates.append(coord_copy)

###############################################################################
# Create ``DpfPath`` object
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a ``DpfPath`` object.

path = post.create_path_on_coordinates(coordinates=coordinates)

###############################################################################
# Request result on path
# ~~~~~~~~~~~~~~~~~~~~~~
# Request the result on this path.

stress = solution.stress(path=path)

###############################################################################
# Plot result
# -----------
# Plot the result.

stress_eqv = stress.von_mises
stress_eqv.plot_contour()
