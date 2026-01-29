"""
.. _ref_cyclic_mesh_skin_example:

Reduce Cyclic Model Size by using Mesh Skin for Result and Mesh extraction
==========================================================================
This example displays post-processing on a mesh skin for a cyclic modal analysis.
The skin mesh is rebuilt with new surface elements connecting the nodes on the external skin
of the solid mesh. These surface elements types are chosen with respect to the solid elements
facets having all their nodes on the skin.
This feature is available for all types of Mechanical simulation supporting cyclic
(or cyclic multistage) and allows you to reduce the size
of the mesh and of the extracted data to improve processing performance. Since larger stress and
strains are usually located on the skin of a model, computing results on the skin gives
equivalent maximum values in most cases. Post-processing of elemental or elemental nodal results
requires element solid to skin mapping to get from a solid element result to a facet
result. The connectivity of the new surface elements built on the skin being different from
the connectivity of the solid elements, small differences can be found after result averaging.
To plot cyclic expanded results, the skin mesh is expanded.
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
# Extract displacement data on the skin.

displacement_skin = simulation.displacement(skin=True)
displacement_skin.plot()


###############################################################################
# Extract stress/strain data
# --------------------------
# Extract stress or elastic strain data over the entire mesh and on the skin.
# Averaging, and invariants computation are done through a solid to skin connectivity mapping.

elemental_stress_skin = simulation.stress_principal_elemental(components=[1], skin=True)
elemental_stress_skin.plot()

elastic_strain_eqv_skin = simulation.elastic_strain_eqv_von_mises_nodal(skin=True)
elastic_strain_eqv_skin.plot()

###############################################################################
# Get stress results on the skin of the first sector with a cyclic phase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stress_eqv_cyc_phase = simulation.stress_eqv_von_mises_nodal(
    set_ids=[5],
    expand_cyclic=[1],
    phase_angle_cyclic=45.0,
    skin=True,
)
stress_eqv_cyc_phase.plot()
