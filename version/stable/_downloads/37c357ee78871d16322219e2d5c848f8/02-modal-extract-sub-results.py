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
.. _ref_modal_sub_results_example:

Extract components of results (modal simulation)
================================================
This example processes a modal simulation to extract subcomponents
of results like displacement and elastic strain.
"""

###############################################################################
# Perform required imports
# ------------------------
# This example uses a supplied file that you can get by importing the DPF ``examples``
# package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.download_all_kinds_of_complexity_modal()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.ModalMechanicalSimulation(example_path)
# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Extract X displacement over a list modes
# ----------------------------------------
# To help pick the right modes, printing the time frequency support.

print(simulation.time_freq_support)

# To get X displacements on the first two modes
x_displacement = simulation.displacement(modes=[1, 2], components=["X"])
# equivalent to
x_displacement = simulation.displacement(set_ids=[1, 2], components=["X"])
print(x_displacement)

x_displacement.plot(set_id=1)

###############################################################################
# Extract XX and XY elastic strain over a list modes
# --------------------------------------------------
# Get X displacements on the first two modes.
XX_XY_elastic_strain = simulation.elastic_strain_nodal(
    modes=[3], components=["XX", "XY"]
)
print(XX_XY_elastic_strain)
