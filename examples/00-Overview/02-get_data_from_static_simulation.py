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
import ansys.dpf.post as dpf
from ansys.dpf.post import examples

simulation = dpf.load_simulation(examples.static_rst)
print(simulation)

###############################################################################
# Get and plot displacements
# --------------------------
displacement = simulation.displacement()

###############################################################################
# Print information
print(displacement._fc)

###############################################################################
# Plot displacements
displacement._fc[0].plot()

###############################################################################
# Get and plot stresses
# ---------------------
# Request "XY" stress component averaged on nodes
stress = simulation.stress_nodal(component_ids="XY")

###############################################################################
# Print information
print(stress._fc)

###############################################################################
# Plot available stresses.
stress._fc[0].plot()

###############################################################################
# Get stresses at only 5 nodes
# ------------------------------
# Request stress only at the first 5 nodes using their IDs.
stress_nodes = simulation.stress_nodal(node_ids=range(1, 6))

###############################################################################
# Print information
print(stress_nodes._fc)

###############################################################################
# Plot stresses
stress_nodes._fc[0].plot()

###############################################################################
# Get stresses in a named selection
# ---------------------------------
# Get the name of the first named selection in the simulation
ns = simulation.named_selections[0]
# Request nodal stresses for this named selection
stress_named_sel = simulation.stress_nodal(named_selections=ns)

###############################################################################
# Print information
print(stress_named_sel._fc)

###############################################################################
# Plot stresses
stress_named_sel._fc[0].plot()

###############################################################################
# Get stresses in a few elements
# ------------------------------
# Request stress only for a few elements selected by their ID
stress_elements = simulation.stress_nodal(element_ids=[1, 2, 3])

###############################################################################
# Print information
print(stress_elements._fc)

###############################################################################
# Plot stresses
stress_elements._fc[0].plot()

###############################################################################
# Get elemental stress and raw stresses
# -------------------------------------
# Request elemental stresses and print information
stress_elements = simulation.stress_elemental()
print(stress_elements._fc)

###############################################################################
# Request raw stresses ("ElementalNodal") and print information
stress_raw = simulation.stress()
print(stress_raw._fc)
