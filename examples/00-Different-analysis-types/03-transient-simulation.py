"""
.. _ref_transient_example:

Transient Simulation with Animation
===================================
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
displacement.animate(deform=True)

# get the available time set ids in the simulation
print(simulation.set_ids)

# extract displacement on given time steps or select the times steps from the already evaluated
# displacement DataFrame
displacement = simulation.displacement(set_ids=simulation.set_ids[5:])
displacement = displacement.select(set_id=simulation.set_ids[5:])
print(displacement)

###############################################################################
# Extract strain at all times or on a selection
# ---------------------------------------------------
strain = simulation.elastic_strain_nodal(all_sets=True)
print(strain)

strain = simulation.elastic_strain_nodal(set_ids=simulation.set_ids[10:])
print(strain)


###############################################################################
# Animate strain eqv over all times
# ---------------------------------

strain_eqv = simulation.elastic_strain_eqv_von_mises_nodal(all_sets=True)
strain_eqv.animate()
