"""
.. _ref_multi-stage_cyclic_results_example:

Extract multi-stage cyclic results
==================================
In this script, a multi-stage modal analysis with cyclic symmetry is processed to show
how to expand the mesh and results.
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

example_path = examples.download_multi_stage_cyclic_result()
simulation = post.ModalMechanicalSimulation(example_path)

#############################################################################
# Extract expanded displacement norm
# ----------------------------------

displacement_norm = simulation.displacement(
    norm=True,
    expand_cyclic=True,
)
print(displacement_norm)
displacement_norm.plot()

#############################################################################
# Extract equivalent von Mises nodal stress without expansion
# -----------------------------------------------------------

stress_vm_sector_0_both_stages = simulation.stress_eqv_von_mises_nodal(
    expand_cyclic=False,
)
print(stress_vm_sector_0_both_stages)
stress_vm_sector_0_both_stages.plot()

#############################################################################
# Extract equivalent von Mises nodal stress expanded on the first four sectors of the first stage
# -----------------------------------------------------------------------------------------------

stress_vm_sectors_0_1_2_3_first_stage = simulation.stress_eqv_von_mises_nodal(
    expand_cyclic=[0, 1, 2, 3],
)
print(stress_vm_sectors_0_1_2_3_first_stage)
stress_vm_sectors_0_1_2_3_first_stage.plot()

#############################################################################
# Extract equivalent von Mises nodal stress expanded on the first two sectors of both stages
# ------------------------------------------------------------------------------------------

stress_vm_sectors_0_1_both_stages = simulation.stress_eqv_von_mises_nodal(
    expand_cyclic=[[0, 1], [0, 1]],
)
print(stress_vm_sectors_0_1_both_stages)
stress_vm_sectors_0_1_both_stages.plot()
