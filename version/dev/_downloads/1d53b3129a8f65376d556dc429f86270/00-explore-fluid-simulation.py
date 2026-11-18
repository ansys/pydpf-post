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
.. _ref_fluid_example_extractions:

Extract fluid simulation results
================================
This example shows how to load a fluid simulation for a heating coil that has been simulated
with Ansys CFX, explore the metadata available, and use the result extraction capabilities.

.. note::
    This example requires DPF 7.0 (2024.1.pre0) or later.
    For more information, see :ref:`compatibility`.

"""
###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get using the ``examples`` module.
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Load fluid analysis result
# --------------------------
fluid_example_files = examples.download_cfx_heating_coil()
simulation = post.FluidSimulation(
    cas=fluid_example_files["cas"], dat=fluid_example_files["dat"]
)
# Print the simulation TO see most of the available metadata
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
# Check available species.
print(simulation.species)

###############################################################################
# Check available phases.
print(simulation.phases)

###############################################################################
# Extract the mesh
# ----------------
###############################################################################
# Extract the full mesh
print(simulation.mesh)
simulation.mesh.plot()
###############################################################################
# Extract the mesh for a given zone by its ID
cell_zone_mesh = simulation.zone_mesh(zone=2)
print(cell_zone_mesh)
cell_zone_mesh.plot()
# or by its name
cell_zone_mesh = simulation.zone_mesh(zone="heater")

face_zone_mesh = simulation.zone_mesh(zone=9)
print(face_zone_mesh)
face_zone_mesh.plot()
face_zone_mesh = simulation.zone_mesh(zone="outflow")

###############################################################################
# Extract a result
# ----------------
###############################################################################
# Check the metadata on available results
print(simulation.result_info)

###############################################################################
# Extract a result
# ----------------
# Print information about a specific available result.
print(simulation.result_info["temperature"])

# Printing information using an index.
print(simulation.result_info[12])

###############################################################################
# Extract temperature data
# -------------------------
temperature = simulation.temperature()
print(temperature)

# The resulting dataframe shows data for two different phases.

###############################################################################
# Select data for phase 2 (Water) only
print(temperature.select(phase=[2]))

###############################################################################
# To directly extract the temperature data for a given phase,
# pass the ``temperature()`` method a ``phases`` argument.
# This argument must be given a list of phase unique IDs, which appear
# in the dataframe in the phase label column between parentheses,
# or as listed under the ``Available qualifier labels`` section of the
# metadata on the result. You can also directly use the phase name.
water_temperature = simulation.temperature(phases=["Water at 25 C"])
# equivalent to
# water_temperature = simulation.temperature(phases=[2])
print(water_temperature)
# # The resulting dataframe now only stores data for the Water phase.

###############################################################################
# To extract a result on given zones use the ``zone_ids`` argument
# or the ``qualifiers`` dictionary argument with the ``zone`` key.
# This code gets and plots the temperature on all face zones.
face_temperature = simulation.temperature(zone_ids=list(simulation.face_zones.keys()))
face_temperature.plot()
