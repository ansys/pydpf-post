"""
.. _ref_fluid_example:

Fluid Simulation
================
This example shows how to load a fluid simulation and extract results like velocity and pressure.
It then shows extraction and selection of results by zone, phase or species.
"""
###############################################################################
# Use gRPC protocol on Linux
# --------------------------
# The CFF plugin is currently prone to errors when used InProcess on Linux,
# hence a gRPC server configuration is chosen when running on a Unix system.
import platform

import ansys.dpf.core as dpf

###############################################################################
# Perform required imports
# ------------------------
from ansys.dpf import post
from ansys.dpf.post import examples

if "Linux" in platform.system():
    dpf.SERVER_CONFIGURATION = dpf.server_factory.AvailableServerConfigs.GrpcServer

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
# Explore the available structure
# -------------------------------
# Check the available cell and face zones
# print(simulation.zones)

###############################################################################
# Check the available species
print(simulation.species)

###############################################################################
# Check the available phases
print(simulation.phases)


###############################################################################
# Extract a result
# ----------------
###############################################################################
# Check the metadata on available results
print(simulation.result_info)

###############################################################################
# You can also get a direct list of ``AvailableResult`` objects
print(simulation.results)
# Print a specific one to get more information on available qualifiers (zones, phases and so on)
print(simulation.results[0])

###############################################################################
# Extract this available result as a Dataframe
print(simulation.enthalpy())
# Not specifying any qualifier returns a unique column of data

###############################################################################
# Available qualifiers for this result can be used in the extraction request
# to filter/separate data.
# Here we want a different column for each available zone for this result
print(simulation.enthalpy(zone_ids=[13, 28]))
# Here only the data for the fluid-stator zone is extracted
print(simulation.enthalpy(zone_ids=[28]))
# The same logic can be applied to any available qualifier found in the ``AvailableResult``
# description under 'Available qualifier labels'.

###############################################################################
# The result extraction request can also contain selection arguments
# The enthalpy result being defined on cells, you can request data for specific
# cells using their IDs:
print(simulation.enthalpy(cell_ids=[1, 2]))

###############################################################################
# For selection and manipulation of the Dataframe,
# please refer to example "Create and manipulate a DPF Dataframe"
