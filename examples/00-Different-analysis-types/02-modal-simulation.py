"""
.. _ref_modal_sub_results_example:

Modal Simulation
================
Simple post processing operations like viewing the different mode shapes is diplayed
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
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
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

for set_id in simulation.set_ids:
    displacement_norm.plot(set_id=set_id)


###############################################################################
# Extract a selection of mode shapes and view them one by one
# -----------------------------------------------------------

modes = [1, 2, 3]

displacement_norm = simulation.displacement(modes=modes, norm=True)
print(displacement_norm)

for set_id in modes:
    displacement_norm.plot(set_id=set_id)
