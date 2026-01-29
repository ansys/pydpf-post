"""
.. _ref_export_data_example:

Export data contained in a DataFrame to a new format
====================================================
In this script a static simulation is used as an example to show how to
export data from a DataFrame to another format.
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

example_path = examples.download_crankshaft()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract elemental nodal stress
# ------------------------------

stress = simulation.stress_nodal(all_sets=True)
print(stress)

###############################################################################
# Export to a numpy.ndarray
# -------------------------

# To export to a numpy.ndarray, the DataFrame must contain data only for a single combination
# of values for column labels.

# select one step as our dataframe contains data for several
stress_1 = stress.select(set_ids=1)
print(stress_1)

# export it as a numpy.ndarray
stress_1_array = stress_1.array
print(stress_1_array)
