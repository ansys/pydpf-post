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
.. _ref_data_data_frame_example:

Explore the data of a result with the DataFrame - Harmonic Simulation
=====================================================================
In this script a harmonic simulation is used as an example to show how to
interact with the post DataFrame: object returned by each result.
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

###############################################################################
# Extract displacement over all sets as an example
# ------------------------------------------------

displacement = simulation.displacement(all_sets=True)
print(displacement)
type(displacement)

###############################################################################
# Loop over all columns and rows to understand the DataFrame
# and get the values for each index.

# columns
for column in displacement.columns:
    print(f'Column with label "{column.name}" and available values {column.values}.')

# rows
for row in displacement.index:
    print(f'Row with label "{row.name}" and available values {row.values}.')


###############################################################################
# Make selections in this DataFrame
# ---------------------------------
# All the labels and values displayed above can be used to select sub parts of the
# DataFrame.

all_real_values = displacement.select(complex=0)
print(all_real_values)

all_imaginary_values = displacement.select(complex=1)
print(all_imaginary_values)

sets_values = displacement.select(set_ids=[1, 2])
print(sets_values)

node_values = displacement.select(node_ids=[3548])
print(node_values)

###############################################################################
# Make selections by index in this DataFrame
# ------------------------------------------
# To select values by index for each label, the iselect method can be used.
# The index to ID order follows what is returned by the values on index method used above.

sets_values = displacement.iselect(set_ids=0)
print(sets_values)

node_values = displacement.iselect(node_ids=[0])
print(node_values)

###############################################################################
# Make multi selections in this DataFrame
# ---------------------------------------

real_values_for_one_set_onde_node = displacement.select(
    node_ids=[3548], set_ids=1, complex=0
)
print(real_values_for_one_set_onde_node)

###############################################################################
# Make selections to plot the DataFrame
# -------------------------------------

displacement.plot(set_ids=1, complex=0)
