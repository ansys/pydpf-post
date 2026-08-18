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
.. _ref_fluid_example:

Postprocess a fluid simulation
==============================
This example shows how to load a fluid simulation, explore the model and its available zones,
species, and phases, and then extract a result.

.. note::
    This example requires DPF 7.0 (2024.1.pre0) or later.
    For more information, see :ref:`compatibility`.

"""
###############################################################################
# Perform required imports
# ------------------------
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Load the fluid analysis result
# ------------------------------
fluid_example_files = examples.download_fluent_axial_comp()
simulation = post.FluidSimulation(
    cas=fluid_example_files["cas"], dat=fluid_example_files["dat"]
)
# Printing the simulation will show most of the available metadata
print(simulation)


###############################################################################
# Explore available metadata
# --------------------------
# Check the available cell and face zones.
print(simulation.cell_zones)
print(simulation.face_zones)

###############################################################################
# The mesh metadata is available separately from the mesh
# as accessing the mesh means loading it.
# Use the ``mesh_info`` property to explore the mesh structure to define queries.
print(simulation.mesh_info)

###############################################################################
# Check available species
# ~~~~~~~~~~~~~~~~~~~~~~~
print(simulation.species)

###############################################################################
# Check available phases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
print(simulation.phases)

###############################################################################
# Check metadata on available results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(simulation.result_info)

###############################################################################
# Extract a result
# ----------------
# Print a specific result to get more information on available qualifiers
# (such as zones and phases).
print(simulation.result_info["enthalpy"])
# Or use an index
# print(simulation.result_info[0])

###############################################################################
# Extract this result as a dataframe
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
enthalpy = simulation.enthalpy()
print(enthalpy)
# Not specifying any qualifier returns a unique column of data

###############################################################################
# Plot dataframe
# ~~~~~~~~~~~~~~
enthalpy.plot()

###############################################################################
# Use qualifiers for result to filter or separate data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# You can use available qualifiers for the result in the extraction request
# to filter or separate data.
#
# This code outputs a different column for each available zone for the result.
print(simulation.enthalpy(zone_ids=[13, 28]))
#
# This code extracts only the data for the fluid-stator zone.
print(simulation.enthalpy(zone_ids=[28]))
#
# You can apply the same logic to any available qualifier found in the ``AvailableResult``
# description under 'Available qualifier labels'.

###############################################################################
# The result extraction request can also contain selection arguments.
# Because the enthalpy result is being defined on cells, you can request data
# for specific cells using their IDs:
print(simulation.enthalpy(cell_ids=[1, 2]))

###############################################################################
# For an example showing how to create and manipulate a DPF dataframe,
# see :ref:`ref_dataframe_example`.
