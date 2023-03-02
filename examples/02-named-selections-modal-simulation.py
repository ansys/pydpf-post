"""
.. _ref_static_example:

Extract results on named selections - Modal Simulation
=======================================================
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

example_path = examples.download_modal_frame()
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Get the available named selections
# ----------------------------------

print(simulation.named_selections)

###############################################################################
# Extract displacements on named selections
# -----------------------------------------

bar1_tot_displacement = simulation.displacement(named_selections=['BAR_1'], norm=True)
print(bar1_tot_displacement)
bar1_tot_displacement.plot()

bar2_tot_displacement = simulation.displacement(named_selections=['BAR_2'], norm=True)
print(bar2_tot_displacement)
bar2_tot_displacement.plot()

# both
tot_displacement = simulation.displacement(named_selections=['BAR_1', 'BAR_2'], norm=True)
print(tot_displacement)
tot_displacement.plot()

###############################################################################
# Extract stress and averaged stress on named selections
# ------------------------------------------------------
eqv_stress = simulation.stress_eqv_von_mises_nodal(named_selections=['_FIXEDSU'])
print(eqv_stress)

# without selection
elemental_stress = simulation.stress_elemental(named_selections=['BAR_1'])
print(elemental_stress)
elemental_stress.plot()
