"""
.. _ref_modal_sub_results_example:

Extract components of results - Modal Simulation
================================================
In this script, a modal simulation is processed to extract sub-components
of results like elastic strain and displacement.
"""

###############################################################################
# Perform required imports
# ------------------------
# This example uses a supplied file that you can
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

example_path = examples.download_all_kinds_of_complexity_modal()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.ModalMechanicalSimulation(example_path)
# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Extract X displacement over a list modes
# ----------------------------------------
# Printing the time freq support can help pick the right modes

print(simulation.time_freq_support)

# To get X displacements on the first 2 modes
x_displacement = simulation.displacement(modes=[1, 2], components=["X"])
# equivalent to
x_displacement = simulation.displacement(set_ids=[1, 2], components=["X"])
print(x_displacement)

x_displacement.plot(set_id=1)

###############################################################################
# Extract XX and XY elastic strain over a list modes
# --------------------------------------------------
# To get X displacements on the first 2 modes
XX_XY_elastic_strain = simulation.elastic_strain_nodal(
    modes=[3], components=["XX", "XY"]
)
print(XX_XY_elastic_strain)
