# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
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
.. _ref_basics:

Get and use a ``Result`` object
===============================
This example shows how to use the legacy PyDPF-Post API to get and use a ``Result`` object.
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
# Get the ``Solution`` object and instantiate with the path for the result
# file.

example_path = examples.download_all_kinds_of_complexity()
solution = post.load_solution(example_path)

###############################################################################
# Get ``Result`` objects
# ----------------------

###############################################################################
# Get displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~

displacement_result = solution.displacement()
displacement = displacement_result.vector

###############################################################################
# Get information on displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

displacement.num_fields
disp_data = displacement.get_data_at_field(0)
len(disp_data)

disp_data[1]

displacement.max_data
displacement.get_max_data_at_field(0)

displacement.min_data

###############################################################################
# Get tensor stress result
# ~~~~~~~~~~~~~~~~~~~~~~~~
# You can get the nodal or elemental location for a tesnsor stress result.
# The default is the nodal location.

el_stress_result = solution.stress(location=post.locations.elemental)
nod_stress_result = solution.stress(location=post.locations.nodal)

###############################################################################
# Get information on tensor result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

el_stress = el_stress_result.tensor
nod_stress = nod_stress_result.tensor

el_field = el_stress[0]
el_field.location

nod_field = nod_stress[0]
nod_field.location

el_stress.get_max_data_at_field(0)
