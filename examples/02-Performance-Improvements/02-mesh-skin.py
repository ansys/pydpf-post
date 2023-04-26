"""
.. _ref_mesh_skin_example:

Reduce Model Size by using Mesh Skin for Result and Mesh extraction
===================================================================
This example displays post-processing on a mesh skin for a static analysis.
The skin mesh is rebuilt with new surface elements connecting the nodes on the external skin
of the solid mesh. These surface elements types are chosen with respect to the solid elements
facets having all their nodes on the skin.
This feature is available for all types of Mechanical simulation and allows you to reduce the size
of the mesh and of the extracted data to improve processing performance. Since larger stress and
strains are usually located on the skin of a model, computing results on the skin gives
equivalent maximum values in most cases. Post-processing of elemental or elemental nodal results
requires element solid to skin mapping to get from a solid element result to a facet
result. The connectivity of the new surface elements built on the skin being different from
the connectivity of the solid elements, small differences can be found after result averaging.
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

example_path = examples.download_crankshaft()
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Extract displacement data
# -------------------------
# Extract displacement data on the skin.

displacement_skin = simulation.displacement(skin=True)
displacement_skin.plot()

print(f"number of nodes with skin=True: {len(displacement_skin.index.mesh_index)}")
print(f"number of nodes with skin=False: {len(simulation.mesh.node_ids)}")

###############################################################################
# Extract stress/strain data
# --------------------------
# Extract stress or elastic strain data and on the skin.
# Averaging, and invariants computation are done through a solid to skin connectivity mapping.

elemental_stress_skin = simulation.stress_principal_elemental(components=[1], skin=True)
elemental_stress_skin.plot()

print(
    f"number of elements with skin=True: {len(elemental_stress_skin.index.mesh_index)}"
)
print(f"number of elements with skin=False: {len(simulation.mesh.element_ids)}")


elastic_strain_eqv_skin = simulation.elastic_strain_eqv_von_mises_nodal(skin=True)
elastic_strain_eqv_skin.plot()

###############################################################################
# Extract the external layer on a selection of elements
# -----------------------------------------------------

all_elements = simulation.mesh.element_ids
elements = []
for i in range(0, int(all_elements.size / 2)):
    elements.append(all_elements[i])
elemental_stress_skin = simulation.stress_principal_elemental(skin=elements)
elemental_stress_skin.plot()
