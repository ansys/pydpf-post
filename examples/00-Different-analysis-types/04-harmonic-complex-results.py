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
.. _ref_harmonic_example:

Harmonic Simulation
===================
In this script harmonic simulation is processed and complex results are used.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
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

example_path = examples.download_harmonic_clamped_pipe()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.HarmonicMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available

print(simulation)


###############################################################################
# Extract displacement over a list of frequencies sets
# ----------------------------------------------------
# Printing the time freq support can help pick the right frequencies

print(simulation.time_freq_support)

displacement = simulation.displacement(set_ids=[1, 2])
print(displacement)

subdisp = displacement.select(complex=0, set_ids=1)
print(subdisp)
subdisp.plot(title="U tot real")

subdisp = displacement.select(complex=1, set_ids=1)
print(subdisp)
subdisp.plot(title="U tot imaginary")

subdisp = displacement.select(complex=0, set_ids=2)
print(subdisp)
subdisp.plot(title="U tot real")

###############################################################################
# Extract stress eqv over a list of frequencies sets
# --------------------------------------------------

stress_eqv = simulation.stress_eqv_von_mises_nodal(set_ids=[1, 2])
print(stress_eqv)

sub_eqv = stress_eqv.select(complex=0, set_ids=1)
print(sub_eqv)
sub_eqv.plot(title="S_eqv real")

sub_eqv = stress_eqv.select(complex=1, set_ids=1)
print(sub_eqv)
sub_eqv.plot(title="S_eqv imaginary")

sub_eqv = stress_eqv.select(complex=0, set_ids=2)
print(sub_eqv)
sub_eqv.plot(title="S_eqv real")
