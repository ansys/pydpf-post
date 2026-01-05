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
.. _ref_modal_example:

Postprocess a modal mechanical simulation
=========================================
This example shows how to postprocess a modal mechanical simulation to view different
mode shapes.
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

###############################################################################
# View the frequency domain
# -------------------------
# To help pick the right modes, printing the time frequency support.

print(simulation.time_freq_support)

# set_ids returns the unique identifiers for the modes
print(simulation.set_ids)

###############################################################################
# Extract and view all mode shapes
# --------------------------------
# Extract all mode shapes and view them one by one.

displacement_norm = simulation.displacement(all_sets=True, norm=True)
print(displacement_norm)

for set_id in simulation.set_ids:
    displacement_norm.plot(set_ids=set_id)


###############################################################################
# Extract and view a selection of mode shapes
# -------------------------------------------
# Extract and view a selection of mode shapes and view them one by one.

modes = [1, 2, 3]

displacement_norm = simulation.displacement(modes=modes, norm=True)
print(displacement_norm)

for set_id in modes:
    displacement_norm.plot(set_ids=set_id)
