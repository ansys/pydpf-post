"""
.. _ref_dataframe_example:

Create and manipulate a DPF Dataframe
=====================================
In this script a Dataframe is generated by extracting a result from a static simulation.
It then showcases different Dataframe viewing and manipulation possibilities.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get the ``Simulation`` object
# -----------------------------
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
# Get a ``Dataframe`` object
# --------------------------
# Extract a result as a Dataframe
displacement_dataframe = simulation.displacement(all_sets=True)

# The Dataframe is displayed as a table, with row and column labels to identify the data.
print(displacement_dataframe)

###############################################################################
# Explore ``Index`` objects
# -------------------------
# The data labels are each defined by an Index object or one of its specialized subtypes.

# The Dataframe's column labels are defined in Dataframe.columns.
print(displacement_dataframe.columns)

###############################################################################
# A ``ResultIndex`` index defines the result stored in the Dataframe.
print(displacement_dataframe.columns[0])
# print(displacement_dataframe.columns.results_index)  # equivalent
###############################################################################
# You can check values available for an Index
print(displacement_dataframe.columns[0].values)

###############################################################################
# A ``SetIndex`` index defines the available set IDs available.
# A set ID is a unique identifier associated to each time-step, step and sub-step, or frequency
# available in a simulation.
# As shown next, an Index has a name and a list of values of a given type.
print(displacement_dataframe.columns[1])
print(displacement_dataframe.columns[1].values)

###############################################################################
# The Dataframe's row labels are defined in Dataframe.index.
print(displacement_dataframe.index)

###############################################################################
# A ``MeshIndex`` defines the mesh entities for which data is available.
# It can store node IDs, element IDs, or face IDs.
print(displacement_dataframe.index[0])
# print(displacement_dataframe.index.mesh_index)  # equivalent
###############################################################################
# Since the list of possible values can be long and querying it can be costly,
# the list of available values may not be determined unless explicitly asked.
print(displacement_dataframe.index[0].values)
###############################################################################
# The ``MeshIndex`` will then be updated to display the actual number of entities available.
print(displacement_dataframe.index[0])
# IMPORTANT: Note that the mesh entity IDs ordered based on the internal data storage structure,
# they are not by ascending order by default!

###############################################################################
# A ``CompIndex`` defines the result components for which data is available.
print(displacement_dataframe.index[1])
print(displacement_dataframe.index[1].values)

###############################################################################
# Change the Dataframe print
# --------------------------
# Options exist to configure the way a Dataframe is displayed.
# You can change the number of data rows displayed with:
displacement_dataframe.display_max_rows = 9
print(displacement_dataframe)
###############################################################################
# Or the number of data columns displayed with:
displacement_dataframe.display_max_columns = 2
print(displacement_dataframe)
# Notice that the ``...`` symbols specify that the Dataframe is truncated in that direction.

###############################################################################
# The special case of ElementalNodal results
# ------------------------------------------
# When dealing with results located on each node of each element (aka. ElementalNodal),
# an ``ElementNodeIndex`` index is added at creation to specify the node number in the
# element's connectivity.
stress = simulation.stress()
print(stress)
print(stress.columns[2])

###############################################################################
# Data selection
# --------------
# To select specific columns or rows, use the index names as arguments for the ``DataFrame.select``
# method, taking lists of values:
disp_X_1 = displacement_dataframe.select(
    set_ids=[1], node_ids=[4872, 9005], components=["X"]
)
print(disp_X_1)

###############################################################################
# You can also select along an index using a zero-based position with ``Dataframe.iselect``:
disp_Y_9005_3 = displacement_dataframe.iselect(
    set_ids=[2], node_ids=[1], components=[1]
)
print(disp_Y_9005_3)

###############################################################################
# Extract data
# ------------
# Once the Dataframe contains the specific data you require, extract it as an array with:
print(disp_X_1.array)
# IMPORTANT: Note that for the extraction of the Dataframe's data as an array to make sense,
# you must first filter the columns label values to a unique combination of values.
# The exception is for ElementalNodal data, which is returned as a 2D array.
print(stress.array.shape)

###############################################################################
# Plot a Dataframe
# ------------------
displacement_dataframe.plot()

###############################################################################
# Animate a transient Dataframe
# -----------------------------
displacement_dataframe.animate()
