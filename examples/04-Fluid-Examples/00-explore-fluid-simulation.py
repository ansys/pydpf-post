"""
.. _ref_fluid_example_extractions:

Fluid Simulation Result Extraction
==================================
This example shows how to load a fluid simulation (here a heating coil simulated with CFX),
explore the metadata available, and use the result extraction capabilities.

.. note::
    This example requires DPF 7.0 (2024.1.pre0) or above.
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
fluid_example_files = examples.download_cfx_heating_coil()
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
# Print information about a specific available result
# for result in simulation.result_info:
#     print(result)
print(simulation.result_info["enthalpy"])
# Or use an index
# print(simulation.result_info[12])

###############################################################################
# Extract the temperature data
# ----------------------------
temperature = simulation.enthalpy(zone_ids=[1, 9])
print(temperature)
print(temperature.select(zone=[9]))
# # The dataframe obtained shows data for two different phases
# print(temperature.select(phase=[2]))
# ###############################################################################
# # To directly extract the temperature data for only one phase,
# # pass the 'temperature' method a 'phases' argument.
# # This argument must be given a list of phase unique identifiers, which appear
# # in the dataframe in the phase label column between parentheses,
# # or as listed
# # under the 'Available qualifier labels' section of the metadata on the result
# # water_temperature = simulation.temperature(phases=["Copper"])
# water_temperature = simulation.enthalpy(
#     phases=[2]
# )
# print(water_temperature)
# # The dataframe obtained now only stores data for the water phase.
