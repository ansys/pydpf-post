"""
.. _ref_get_data_from_transient_simulation:

Get data from transient simulation
==================================
This example show how to request data from a previously stored transient simulation.
The available results can be listed to see what results can be retrieved.
"""

###############################################################################
# Imports and loading simulation
# ------------------------------
from ansys.dpf.post import examples, load_simulation

simulation = load_simulation(examples.download_transient_result())
print(simulation)

###############################################################################
# Get displacements for all timesteps
# -----------------------------------
disp = simulation.displacement()

###############################################################################
# Print displacement information
print(disp)

###############################################################################
# Get information for the last timestep displacement field
print(disp[-1])

###############################################################################
# Get data for the last timestep displacement field
print(disp[-1].data)

###############################################################################
# Plot data for the last timestep displacement field
disp[-1].plot()

###############################################################################
# Get displacements for one timestep
# -----------------------------------
# Request displacements only for the selectet time steps. Note that ``steps`` argument
# allow you to select a certain time-step by its index (integers from 1 to n_time_steps).
# This option is more efficient if only a few time steps are required.
disp_times = simulation.displacement(steps=[2, 22])

###############################################################################
# Print displacement information
print(disp_times)

###############################################################################
# Get nodal stresses for a few nodes
# ----------------------------------
nodal_stress = simulation.nodal_stress(nodes=range(1000, 1400))

###############################################################################
# Print information. Note that only 58 nodes where found with IDs between
# 1000 and 1400.
print(nodal_stress)

###############################################################################
# Select the first time step and plot nodal stress in the requested nodes.
nodal_stress[0].plot()

###############################################################################
# Select the last time step and plot nodal stress in the requested nodes.
nodal_stress[-1].plot()
