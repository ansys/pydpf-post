"""
.. _ref_cyclic_mesh_external_layer_example:

Reduce Cyclic Model Size by using Mesh External Layer for Result and Mesh extraction
====================================================================================
This example displays post-processing on a mesh external layer for a cyclic modal analysis.
The external layer is the layer of solid elements with at least one facet facing the outside of
the geometry.
This feature is available for all types of Mechanical simulation supporting cyclic
(or cyclic multistage) and allows you to reduce the size
of the mesh and of the extracted data to improve processing performance.
Since larger stress and strains are usually located on the skin of a model,
computing results on the external layer gives equivalent maximum values in most cases.
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

example_path = examples.download_modal_cyclic()
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract displacement data
# -------------------------
# Extract displacement data on the external layer.

displacement_ext = simulation.displacement(external_layer=True)
print(displacement_ext)
displacement_ext.plot()


###############################################################################
# Select sectors for expansion
# ----------------------------
# Extract displacement data on the external layer
# for a cyclic expansion on a selected number of sectors.

displacement_ext = simulation.displacement(external_layer=True, expand_cyclic=[1, 2, 3])
displacement_ext.plot()

###############################################################################
# Extract stress/strain data
# --------------------------
# Extract stress or elastic strain data on the external layer.
# Averaging, and invariants computation can easily be done on the external layer since the
# connectivity of the kept elements remains unchanged.

elemental_stress_ext = simulation.stress_principal_elemental(
    components=[1], external_layer=True
)
elemental_stress_ext.plot()

elastic_strain_eqv_ext = simulation.elastic_strain_eqv_von_mises_nodal(
    external_layer=True
)
elastic_strain_eqv_ext.plot()

##################################################################################
# Get stress results on the external layer of the first sector with a cyclic phase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stress_eqv_cyc_phase = simulation.stress_eqv_von_mises_nodal(
    set_ids=[5],
    expand_cyclic=[1],
    phase_angle_cyclic=45.0,
    external_layer=True,
)
stress_eqv_cyc_phase.plot()
