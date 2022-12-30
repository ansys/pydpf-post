"""
.. _ref_get_data_from_static_simulation:

Get data from static simulation
===============================
This example shows how to request data from a previously stored static simulation.
The available results can be listed to see what results can be retrieved.
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

###############################################################################
# Print information
print(disp)

###############################################################################
# Plot displacements
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

###############################################################################
# Get stresses in a named selection
# ---------------------------------
# Request stress in named selection "_CM82"
stress_named_sel = simulation.nodal_stress(named_selection="_CM82")

###############################################################################
# Print information
print(stress_named_sel)

###############################################################################
# Plot stresses
stress_named_sel[0].plot()

###############################################################################
# Get stresses in a few elements
# ------------------------------
# Request stress only for a few elements selected by their ID
stress_elements = simulation.nodal_stress(elements=[301, 302, 303])

###############################################################################
# Print information
print(stress_elements)

###############################################################################
# Plot stresses
stress_elements[0].plot()

###############################################################################
# Get elemental stress and raw stresses
# -------------------------------------
# Request elemental stresses and print information
stress_elements = simulation.elemental_stress()
print(stress_elements)

###############################################################################
# Request raw stresses ("ElementalNodal") and print information
stress_raw = simulation.raw_stress()
print(stress_raw)
