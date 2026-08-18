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
.. _ref_static_example:

Postprocess a static mechanical simulation
==========================================
This example shows how to postprocess a static mechanical simulation to extract results like
displacement and stress. It shows how to selecting subparts of the results by scoping
on specific nodes or elements.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.find_static_rst()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

displacement = simulation.displacement()
print(displacement)


###############################################################################
# Select subparts of displacement
# -------------------------------

# To get X displacements
x_displacement = displacement.select(components="X")
print(x_displacement)


# equivalent to
x_displacement = simulation.displacement(components=["X"])
print(x_displacement)

# plot
x_displacement.plot()

# extract displacement on specific nodes
nodes_displacement = displacement.select(node_ids=[1, 10, 100])
nodes_displacement.plot()

# equivalent to:
nodes_displacement = simulation.displacement(node_ids=[1, 10, 100])
print(nodes_displacement)

###############################################################################
# Compute total displacement (norm)
# ---------------------------------
# Compute the norm of the displacement on a selection of nodes.

nodes_displacement = simulation.displacement(
    node_ids=simulation.mesh.node_ids[10:], norm=True
)
print(nodes_displacement)
nodes_displacement.plot()


###############################################################################
# Extract tensor stresses
# ------------------------
# Extract raw elemental nodal stresses from the result file. Then, apply averaging
# and compute equivalent stresses.
elem_nodal_stress = simulation.stress()
print(elem_nodal_stress)

# Compute nodal stresses from the result file
nodal_stress = simulation.stress_nodal()
print(nodal_stress)

# Compute elemental stresses from the result file
elemental_stress = simulation.stress_elemental()
print(elemental_stress)

# Extract elemental stresses on specific elements
elemental_stress = elemental_stress.select(element_ids=[5, 6, 7])
elemental_stress.plot()

# Compute nodal eqv stresses from the result file
eqv_stress = simulation.stress_eqv_von_mises_nodal()
print(eqv_stress)
eqv_stress.plot()
