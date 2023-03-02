"""
.. _ref_static_example:

Static Simulation
=================
In this script static simulation is processed to extract results like stress, displacement.
Selecting sub parts of the results by scoping on specific nodes, elements is also displayed here.
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

example_path = examples.find_static_rst()
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

displacement = simulation.displacement()
print(displacement)


###############################################################################
# Select sub parts of displacement
# ---------------------------------

# To get X displacements
x_displacement = displacement.select(comp="X")
print(x_displacement)


# equivalent to 
x_displacement = simulation.displacement(components=["X"])
print(x_displacement)

# plot
x_displacement.plot()

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
