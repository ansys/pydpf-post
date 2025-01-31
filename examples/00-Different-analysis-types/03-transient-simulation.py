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
.. _ref_transient_example:

Postprocess a transient mechanical simulation with animation
============================================================
This example shows how to postprocess a transient mechanical simulation with animation
to extract results like displacement, stress, and strain. It also shows how
to extract data for chosen time steps and animate the strain equivalence over
all times.
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

example_path = examples.find_msup_transient()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.TransientMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Extract displacement
# --------------------
# You can extract displacement at all times or on a selection of time steps.

# query the displacement vectorial field for all times
displacement = simulation.displacement(all_sets=True)
print(displacement)
# animation shows the norm of vectorial fields with several components
displacement.animate(deform=True, title="U")


# get specific components with "components"
x_displacement = simulation.displacement(all_sets=True, components=["X"])
print(x_displacement)
x_displacement.animate(deform=True, title="UX")


# get the norm of a vectorial result with "norm=True"
displacement_norm = simulation.displacement(all_sets=True, norm=True)
print(displacement_norm)
displacement_norm.animate(deform=True, title="U norm")

# get the available time set ids in the simulation
print(simulation.set_ids)

# extract displacement on given time steps or select the times steps from the already evaluated
# displacement DataFrame
displacement = simulation.displacement(set_ids=simulation.set_ids[5:])
displacement = displacement.select(set_ids=simulation.set_ids[5:])
print(displacement)

###############################################################################
# Extract strain
# --------------
# You can extract strain at all times or on a selection of time steps.

strain = simulation.elastic_strain_nodal(all_sets=True)
print(strain)

strain = simulation.elastic_strain_nodal(set_ids=simulation.set_ids[10:])
print(strain)


###############################################################################
# Animate strain eqvivalent over all times
# ----------------------------------------

strain_eqv = simulation.elastic_strain_eqv_von_mises_nodal(all_sets=True)
strain_eqv.animate(deform=True, title="E_eqv")
