"""
.. _ref_get_data_from_simulation:

Get data from simulation
========================
This example show how to request data from a previously stored simulation. The
available results can be listed to see what results can be retrieved.
"""

###############################################################################
# Imports and loading simulation
# ------------------------------
from ansys.dpf.post import examples, load_simulation

simulation = load_simulation(examples.download_all_kinds_of_complexity())

###############################################################################
# Get and plot displacements
# --------------------------
disp = simulation.displacement()

print(disp)

disp[0].plot()

###############################################################################
# Get and plot stresses
# ---------------------
# Request "XY" stress component
stress = simulation.nodal_stress(component="XY")

###############################################################################
# Print information
print(stress)

###############################################################################
# Plot available stresses
stress[0].plot()
stress[1].plot()

###############################################################################
# Get stresses in only 100 nodes
# ------------------------------
# Request stress only in 100 nodes by its ID
stress_nodes = simulation.nodal_stress(nodes=range(400, 500))

###############################################################################
# Print information
print(stress_nodes)

###############################################################################
# Plot stresses
stress_nodes[0].plot()
