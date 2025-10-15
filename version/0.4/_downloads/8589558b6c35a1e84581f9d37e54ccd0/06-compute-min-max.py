"""
.. _ref_compute_statistics_example:

Compute minimum and maximum of a DataFrame
==========================================
In this example, transient mechanical displacement data is used
to show how to compute the min or max of a given DataFrame.
"""

###############################################################################
# Perform required imports
# ------------------------
# This example uses a supplied file that you can
# get using the ``examples`` module.

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
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract displacement data
# -------------------------

displacement = simulation.displacement(all_sets=True)
print(displacement)

###############################################################################
# Compute the maximum displacement for each component at each time-step
# ---------------------------------------------------------------------

# The default axis is the MeshIndex
maximum_over_mesh = displacement.max()
print(maximum_over_mesh)
# is equivalent to
maximum_over_mesh = displacement.max(axis="node_ids")
print(maximum_over_mesh)

# Compute the maximum displacement for each node and component across time
# ------------------------------------------------------------------------
maximum_over_time = displacement.max(axis="set_ids")
print(maximum_over_time)

# Compute the maximum displacement overall
# ----------------------------------------
maximum_overall = maximum_over_time.max()
print(maximum_overall)

###############################################################################
# Compute the minimum displacement for each component at each time-step
# ---------------------------------------------------------------------

# The default axis is the MeshIndex
minimum_over_mesh = displacement.min()
print(minimum_over_mesh)
# is equivalent to
minimum_over_mesh = displacement.min(axis="node_ids")
print(minimum_over_mesh)

# Compute the minimum displacement for each node and component across time
# ------------------------------------------------------------------------
minimum_over_time = displacement.min(axis="set_ids")
print(minimum_over_time)

# Compute the minimum displacement overall
# ----------------------------------------
minimum_overall = minimum_over_time.min()
print(minimum_overall)
