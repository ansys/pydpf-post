"""
.. _ref_data_data_frame_example:

Explore the data of a result with the dataframe (harmonic simulation)
=====================================================================
This example uses a harmonic simulation to show how to interact with the PyDPF-Post
dataframe, which is the object returned by each result.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
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

example_path = examples.download_harmonic_clamped_pipe()
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.HarmonicMechanicalSimulation(example_path)

###############################################################################
# Extract displacement over all sets
# ----------------------------------

displacement = simulation.displacement(all_sets=True)
print(displacement)
type(displacement)

###############################################################################
# Loop over all columns and rows to understand the dataframe
# and get the values for each index.

# columns
for column in displacement.columns:
    print(f'Column with label "{column.name}" and available values {column.values}.')

# rows
for row in displacement.index:
    print(f'Row with label "{row.name}" and available values {row.values}.')


###############################################################################
# Make selections in dataframe
# ----------------------------
# Use the labels and values shown by the preceding code to select subparts of the
# dataframe.

all_real_values = displacement.select(complex=0)
print(all_real_values)

all_imaginary_values = displacement.select(complex=1)
print(all_imaginary_values)

sets_values = displacement.select(set_ids=[1, 2])
print(sets_values)

node_values = displacement.select(node_ids=[3548])
print(node_values)

###############################################################################
# Make selections in dataframe by index
# --------------------------------------
# To select values by index for each label, use the ``iselect()`` method.
# The index to ID order follows what is returned by the values on the preceding index method.

sets_values = displacement.iselect(set_ids=0)
print(sets_values)

node_values = displacement.iselect(node_ids=[0])
print(node_values)

###############################################################################
# Make multiple selections in dataframe
# -------------------------------------

real_values_for_one_set_onde_node = displacement.select(
    node_ids=[3548], set_ids=1, complex=0
)
print(real_values_for_one_set_onde_node)

###############################################################################
# Make selections to plot dataframe
# ---------------------------------

displacement.plot(set_ids=1, complex=0)
