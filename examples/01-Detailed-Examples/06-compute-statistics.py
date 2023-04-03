"""
.. _ref_compute_statistics_example:

Compute statistics of a DataFrame
=================================
In this example, transient mechanical displacement data is used
to show how to compute statistical values such as the min or max
for a given DataFrame.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports.
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
# Compute the maximum displacement value for each component at each time-step
# ---------------------------------------------------------------------------

# The default axis is the MeshIndex
maximum = displacement.max()
print(maximum)
# is equivalent to
maximum = displacement.max(axis="node_ids")
print(maximum)

###############################################################################
# Compute the minimum displacement for each node and component across time
# ------------------------------------------------------------------------

minimum = displacement.min(axis="set_ids")
print(minimum)
