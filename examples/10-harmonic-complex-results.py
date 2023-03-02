"""
.. _ref_harmonic_complex_example:

Handling Complex Results - Harmonic Simulation
=============================================
In this script harmonic simulation is processed and complex results are used.
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

example_path = r'D:\ANSYSDev\DPFUnitTestDataFiles\harmonic\RstFile_ENGHRM_WB2_002.rst'
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
simulation = post.HarmonicMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available

print(simulation)


###############################################################################
# Extract displacement over a list of frequencies sets
# ----------------------------------------------------
# Printing the time freq support can help pick the right frequencies

print(simulation.time_freq_support)

displacement = simulation.displacement(set_ids=[1, 2])
print(displacement)

subdisp = displacement.select(complex=0, time=1)
subdisp.plot()
subdisp = displacement.select(complex=1, time=1)
subdisp.plot()
subdisp = displacement.select(complex=1, time=2)
subdisp.plot()