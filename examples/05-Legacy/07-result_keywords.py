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
.. _ref_result_keywords:

Use keywords to get a more precise result
=========================================
This example shows how to use the legacy PyDPF-Post API to use keywords to get
a more precise result when calling a ``Result`` object from a ``Solution`` object.
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
# Get the ``Solution`` object. This example loads a simple file supplied with
# PyDPF-Post.

solution = post.load_solution(examples.multishells_rst)

###############################################################################
# Get keyword list
# ----------------

post.print_available_keywords()


###############################################################################
# Get ``Result`` objects
# ----------------------

###############################################################################
# Get displacement result using scoping
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The default location is nodal.

displacement_result = solution.displacement(
    location=post.locations.nodal, node_scoping=[1, 2, 3]
)
displacement = displacement_result.vector

###############################################################################
# Get information on displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

displacement.get_data_at_field(0)

stress_with_elem_scop_result = solution.stress(
    location=post.locations.elemental_nodal, element_scoping=[1]
)
stress_with_elem_scop = stress_with_elem_scop_result.tensor

stress_with_elem_scop.get_data_at_field(0)

###############################################################################
# Use named selection on displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stress_on_ns_result = solution.stress(
    location=post.locations.elemental_nodal, named_selection="SELECTION"
)
stress_on_ns = stress_on_ns_result.tensor
stress_on_ns.num_fields

len(stress_on_ns[0])

###############################################################################
# Get a subresult
# ~~~~~~~~~~~~~~~

disp_x = displacement_result.x

stress_yz = stress_with_elem_scop_result.yz

stress_principal_1 = stress_on_ns_result.principal_3
stress_principal_1

###############################################################################
# Filter result
# ~~~~~~~~~~~~~
# Filter the result on a time, time scoping, and set.

print(solution.time_freq_support)

stress_on_time_1s_result = solution.stress(time=1.0)
stress_on_time_1s = stress_on_time_1s_result.tensor

displacement_on_set_1_result = solution.displacement(set=1)
displacement_on_set_1 = displacement_on_set_1_result.vector

elastic_strain_with_time_scoping_result = solution.elastic_strain(time_scoping=[1, 3])
elastic_strain_with_time_scoping = elastic_strain_with_time_scoping_result.tensor
elastic_strain_with_time_scoping

###############################################################################
# Make grouping
# ~~~~~~~~~~~~~

displacement_result = solution.displacement(grouping=post.grouping.by_el_shape)
displacement_by_el_shape = displacement_result.vector

###############################################################################
# Filter for MAPDL elements
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter for only MAPDL elements of type solid 186.

stress_result = solution.stress(mapdl_grouping=186)
stress_on_solid_186 = stress_result.tensor
stress_on_solid_186

###############################################################################
# Manipulate result and change its definition
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(stress_on_ns_result)
print(stress_on_ns_result.definition.location)

stress_on_ns_result.definition.location = post.locations.elemental
stress_on_ns_result.definition.time = 1.0
stress_on_ns_elemental = stress_on_ns_result.tensor

print(stress_on_ns_result)

###############################################################################
# Use miscellaneous results
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# You can use the same keywords for miscellaneouse results. For a complex
# result, you can also use the keyword ``phase``, which has a float value.

stress_ratio = solution.misc.elemental_stress_ratio(node_scoping=[1, 32], time=1.0)
print(stress_ratio)
