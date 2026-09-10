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

# noqa: D400
"""
.. _ref_basic_cyclic:

Get base and duplicate sector results for a modal cyclic symmetry model
=======================================================================

This example shows how to extract real and imaginary results from a modal cyclic symmetry model.

"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get using the ``examples`` module.
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.download_modal_cyclic()
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Get base and duplicate sector displacement results
# --------------------------------------------------
# With ``expand_cyclic=False``, the result providers (such as stress and displacement) return
# results for base and duplicate sectors for the cyclic symmetry model.

# Create displacement operator
u_cyc = simulation.displacement(all_sets=True, expand_cyclic=False)

###############################################################################
# The output dataframe print displays the organization of the data.
# The ``base_sector`` label gives access to base sector
# results with ``base_sector=True`` and duplicate sector results with
# ``base_sector=False`` for all modes.
#
# The print also shows that there are no duplicate sectors for the first six modes.
# Indeed, modes with harmonic index 0 have 0.0 displacement, stresses... on
# duplicate sectors.

# print(u_cyc)
print(simulation.time_freq_support)

s_cyc = simulation.stress_eqv_von_mises_nodal(all_sets=True, expand_cyclic=False)
# plot mode 7 base sector (real) result
s_cyc.plot(set_ids=7, base_sector=1)

# plot mode 7 duplicate sector (imaginary) result
s_cyc.plot(set_ids=7, base_sector=0)


###############################################################################
# Get displacement results on first sector with a cyclic phase
# ------------------------------------------------------------
# Get displacemetn results with phi set to different values.

# with phi=0°
u_cyc = simulation.displacement(all_sets=True, expand_cyclic=[1])
u_cyc.plot()

# with phi=90°
u_cyc = simulation.displacement(
    all_sets=True, expand_cyclic=[1], phase_angle_cyclic=90.0
)
u_cyc.plot()

# with phi=45°
u_cyc = simulation.displacement(
    all_sets=True, expand_cyclic=[1], phase_angle_cyclic=45.0
)
u_cyc.plot()

###############################################################################
# Get nodal stress results on first sector with a cyclic phase
# ------------------------------------------------------------
s_cyc = simulation.stress_eqv_von_mises_nodal(
    all_sets=True, expand_cyclic=[1], phase_angle_cyclic=45.0
)
print(s_cyc)
s_cyc.plot()

###############################################################################
# Get elemental nodal stress results on first sector with a cyclic phase
# ----------------------------------------------------------------------
# Elemental nodal is the default result location for stress and strain.
s_cyc = simulation.stress(set_ids=[7], expand_cyclic=[1], phase_angle_cyclic=45.0)
print(s_cyc)

# To average the result for each element
to_elemental = simulation.stress_elemental(
    set_ids=[7], expand_cyclic=[1], phase_angle_cyclic=45.0
)
print(to_elemental)
to_elemental.plot()

###############################################################################
# Get nodal stress results expanded
# ---------------------------------

s_cyc = simulation.stress_eqv_von_mises_nodal(set_ids=[7])
s_cyc.plot()
