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
.. _ref_get_data_from_static_simulation:

Get data from a static simulation
=================================
This example shows how to get data from a previously stored static simulation.
The available results can be listed to see what results can be retrieved.
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
# Get the ``Solution`` object that allows access to the result. The ``Solution``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

simulation = post.load_simulation(examples.static_rst)
print(simulation)

###############################################################################
# Get and plot displacements
# --------------------------
displacement = simulation.displacement()

###############################################################################
# Print information
print(displacement._fc)

###############################################################################
# Plot displacements
displacement._fc[0].plot()

###############################################################################
# Get and plot stresses
# ---------------------
# Get the "XY" stress component averaged on nodes.
stress = simulation.stress_nodal(components="XY")

###############################################################################
# Print information.
print(stress._fc)

###############################################################################
# Plot available stresses.
stress._fc[0].plot()

###############################################################################
# Get stresses at only five nodes
# -------------------------------
# Get stresses at only the first five nodes using their IDs.
stress_nodes = simulation.stress_nodal(node_ids=range(1, 6))

###############################################################################
# Print information.
print(stress_nodes._fc)

###############################################################################
# Plot stresses.
stress_nodes._fc[0].plot()

###############################################################################
# Get stresses in a named selection
# ---------------------------------
# Get the name of the first named selection in the simulation.
ns = simulation.named_selections[0]
# Get nodal stresses for this named selection
stress_named_sel = simulation.stress_nodal(named_selections=ns)

###############################################################################
# Print information.
print(stress_named_sel._fc)

###############################################################################
# Plot stresses.
stress_named_sel._fc[0].plot()

###############################################################################
# Get stresses in a few elements
# ------------------------------
# Get stresses only for a few elements selected by their IDs.
stress_elements = simulation.stress_nodal(element_ids=[1, 2, 3])

###############################################################################
# Print information.
print(stress_elements._fc)

###############################################################################
# Plot stresses.
stress_elements._fc[0].plot()

###############################################################################
# Get elemental streses and raw stresses
# --------------------------------------
# Get elemental stresses and print information.
stress_elements = simulation.stress_elemental()
print(stress_elements._fc)

###############################################################################
# Get raw stresses (ElementalNodal) and print information.
stress_raw = simulation.stress()
print(stress_raw._fc)
