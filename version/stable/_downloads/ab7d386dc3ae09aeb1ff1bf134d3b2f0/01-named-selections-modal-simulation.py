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
.. _ref_ns_modal_example:

Extract results on named selections (modal simulation)
======================================================
This example shows how to process a modal simulation to extract results
like displacement and stress. It selects subparts of the results by scoping
on specific nodes and also shows elements.
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

example_path = examples.download_modal_frame()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Get available named selections
# ------------------------------

print(simulation.named_selections)

###############################################################################
# Extract displacements on named selections
# -----------------------------------------

bar1_tot_displacement = simulation.displacement(named_selections=["BAR_1"], norm=True)
print(bar1_tot_displacement)
bar1_tot_displacement.plot()

bar2_tot_displacement = simulation.displacement(named_selections=["BAR_2"], norm=True)
print(bar2_tot_displacement)
bar2_tot_displacement.plot()

# both
tot_displacement = simulation.displacement(
    named_selections=["BAR_1", "BAR_2"], norm=True
)
print(tot_displacement)
tot_displacement.plot()

###############################################################################
# Extract stress and averaged stress on named selections
# ------------------------------------------------------
eqv_stress = simulation.stress_eqv_von_mises_nodal(named_selections=["_FIXEDSU"])
print(eqv_stress)

# without selection
elemental_stress = simulation.stress_elemental(named_selections=["BAR_1"])
print(elemental_stress)
elemental_stress.plot()
