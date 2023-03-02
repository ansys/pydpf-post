"""
.. _ref_static_example:

Transient Simulation
====================
In this script transient simulation is processed to extract results like
stress, strain, displacement.
Extracting data for chosen time steps and animating is also displayed.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.find_msup_transient()
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
simulation = post.TransientMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Extract displacement at all times or on a selection
# ---------------------------------------------------

displacement = simulation.displacement(all_sets=True)
print(displacement)
displacement.animate(deform=True)


# equivalent to 
x_displacement = simulation.displacement(all_sets=True, components=["X"])
print(x_displacement)

# plot
x_displacement.plot(set_id=1)
x_displacement.animate()

# extract displacement on specific nodes
nodes_displacement = displacement.select(node=[1, 10, 100])
nodes_displacement.plot()

# equivalent to:
nodes_displacement = simulation.displacement(node_ids=[1, 10, 100])
print(nodes_displacement)

###############################################################################
# Compute total displacement (norm)
# ---------------------------------
# Compute the norm of displacement on a selection of nodes
nodes_displacement = simulation.displacement(node_ids=[1, 10, 100], norm=True)
print(nodes_displacement)
nodes_displacement.plot()


###############################################################################
# Extract tensor stress, apply averaging, compute equivalent
# ----------------------------------------------------------
# Extract raw elemental nodal stresses from the rst file
elem_nodal_stress = simulation.stress()
print(elem_nodal_stress)

# Compute nodal stresses from the result file
nodal_stress = simulation.stress_nodal()
print(nodal_stress)

# Compute elemental stresses from the result file
elemental_stress = simulation.stress_elemental()
print(elemental_stress)

# Compute nodal eqv stresses from the result file
eqv_stress = simulation.stress_eqv_von_mises_nodal()
print(eqv_stress)
eqv_stress.plot()
