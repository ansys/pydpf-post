"""
.. _ref_compute_statistics_example:

Compute the minimum and maximum of a dataframe
==============================================
This example uses transient mechanical displacement data
to show how to compute the minimum and maximum of a dataframe.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
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
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract displacement data
# -------------------------

displacement = simulation.displacement(all_sets=True)
print(displacement)

###############################################################################
# Compute maximum displacement for each component at each time-step
# -----------------------------------------------------------------

# The default axis is the MeshIndex
maximum_over_mesh = displacement.max()
print(maximum_over_mesh)
# is equivalent to
maximum_over_mesh = displacement.max(axis="node_ids")
print(maximum_over_mesh)

# Compute maximum displacement for each node and component across time
# --------------------------------------------------------------------
maximum_over_time = displacement.max(axis="set_ids")
print(maximum_over_time)

# Compute maximum displacement overall
# ----------------------------------------
maximum_overall = maximum_over_time.max()
print(maximum_overall)

###############################################################################
# Compute minimum displacement for each component at each time-step
# -----------------------------------------------------------------

# The default axis is the MeshIndex
minimum_over_mesh = displacement.min()
print(minimum_over_mesh)
# is equivalent to
minimum_over_mesh = displacement.min(axis="node_ids")
print(minimum_over_mesh)

# Compute minimum displacement for each node and component across time
# --------------------------------------------------------------------
minimum_over_time = displacement.min(axis="set_ids")
print(minimum_over_time)

# Compute minimum displacement overall
# ------------------------------------
minimum_overall = minimum_over_time.min()
print(minimum_overall)
