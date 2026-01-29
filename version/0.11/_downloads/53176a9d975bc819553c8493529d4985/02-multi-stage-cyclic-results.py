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
.. _ref_multi-stage_cyclic_results_example:

Extract multi-stage cyclic results
==================================
This example uses a multi-stage modal analysis with cyclic symmetry to show
how to expand the mesh and results.
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

example_path = examples.download_multi_stage_cyclic_result()
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

#############################################################################
# Extract expanded displacement normal
# ------------------------------------

displacement_norm = simulation.displacement(
    norm=True,
    expand_cyclic=True,
)
print(displacement_norm)
displacement_norm.plot()

#############################################################################
# Extract equivalent von Mises nodal stress without expansion
# -----------------------------------------------------------

stress_vm_sector_1_both_stages = simulation.stress_eqv_von_mises_nodal(
    expand_cyclic=False,
)
print(stress_vm_sector_1_both_stages)
stress_vm_sector_1_both_stages.plot()

#############################################################################
# Extract equivalent von Mises nodal stress expanded on first four sectors of first stage
# ----------------------------------------------------------------------------------------

stress_vm_sectors_1_2_3_4_first_stage = simulation.stress_eqv_von_mises_nodal(
    expand_cyclic=[1, 2, 3, 4],
)
print(stress_vm_sectors_1_2_3_4_first_stage)
stress_vm_sectors_1_2_3_4_first_stage.plot()

#############################################################################
# Extract equivalent von Mises nodal stress expanded on first two sectors of both stages
# --------------------------------------------------------------------------------------

stress_vm_sectors_1_2_both_stages = simulation.stress_eqv_von_mises_nodal(
    expand_cyclic=[[1, 2], [1, 2]],
)
print(stress_vm_sectors_1_2_both_stages)
stress_vm_sectors_1_2_both_stages.plot()
