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
# This code uses a different column for each available zone for the result.
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
