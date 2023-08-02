"""
.. _ref_modal_example:

Modal Simulation
================
Simple post processing operations like viewing the different mode shapes is displayed
in this example.
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
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.ModalMechanicalSimulation(example_path)

###############################################################################
# View the frequency domain
# -------------------------
# Printing the time freq support can help pick the right modes

print(simulation.time_freq_support)

# set_ids returns the unique identifiers for the modes
print(simulation.set_ids)

###############################################################################
# Extract all mode shapes and view them one by one
# ------------------------------------------------

displacement_norm = simulation.displacement(all_sets=True, norm=True)
print(displacement_norm)

###############################################################################
# Available modes
print(simulation.set_ids)

###############################################################################
# Plot mode 1
displacement_norm.plot(set_ids=1)

###############################################################################
# Plot mode 2
displacement_norm.plot(set_ids=2)

###############################################################################
# Plot mode 3
displacement_norm.plot(set_ids=3)

###############################################################################
# Plot mode 4
displacement_norm.plot(set_ids=4)

###############################################################################
# Plot mode 5
displacement_norm.plot(set_ids=5)

###############################################################################
# Plot mode 6
displacement_norm.plot(set_ids=6)


###############################################################################
# Extract a selection of mode shapes
# ----------------------------------
modes = [1, 2, 3]
displacement_norm = simulation.displacement(modes=modes, norm=True)
print(displacement_norm)
