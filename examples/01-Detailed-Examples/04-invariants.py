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
.. _ref_invariants_example:

Extract stress and strain invariants (static simulation)
========================================================
This example uses a static simulation to show how to extract a tensor's
stress and strain invariants.
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

example_path = examples.download_crankshaft()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract elemental nodal stress and strain
# -----------------------------------------

stress = simulation.stress(all_sets=True)
print(stress)

strain = simulation.elastic_strain(all_sets=True)
print(strain)

###############################################################################
# Compute principal invariant averaged and unaveraged on stress and strain
# ------------------------------------------------------------------------

princ_stress_1 = simulation.stress_principal(components=[1])
print(princ_stress_1)

nodal_princ_stress_2 = simulation.stress_principal_nodal(components=[2])
print(nodal_princ_stress_2)
nodal_princ_stress_2.plot()

nodal_princ_strain_2 = simulation.elastic_strain_principal_nodal(components=[2])
print(nodal_princ_strain_2)
nodal_princ_strain_2.plot()


###############################################################################
# Compute von Mises eqvivalent averaged and unaveraged on stress and strain
# -------------------------------------------------------------------------

stress_eqv = simulation.stress_eqv_von_mises(set_ids=[1])
print(stress_eqv)

nodal_stress_eqv = simulation.stress_eqv_von_mises_nodal(set_ids=[1])
nodal_stress_eqv.plot()

nodal_strain_eqv = simulation.elastic_strain_eqv_von_mises_nodal(set_ids=[1])
nodal_strain_eqv.plot()
