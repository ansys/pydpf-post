"""
.. _ref_mesh_external_layer_example:

Reduce Model Size by using Mesh External Layer for Result and Mesh extraction
=============================================================================
This example displays post-processing on a mesh external layer for a static analysis.
The external layer is the layer of solid elements with at least one facet facing the outside of
the geometry.
This feature is available for all types of Mechanical simulation, and allows you to reduce the size
of the mesh and of the extracted data to improve processing performance.
Since larger stress and strains are usually located on the skin of a model,
computing the results on the external layer provides equivalent maximum values in most cases.
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
# to automatically detect the simulation type, use:
simulation = post.load_simulation(example_path)

# to enable auto-completion, use the equivalent:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract displacement data
# -------------------------
# Extract displacement data on the external layer.

displacement_ext = simulation.displacement(
    external_layer=True
)  # default is external_layer=False
displacement_ext.plot()

print(
    f"number of nodes with external_layer=True: {len(displacement_ext.index.mesh_index)}"
)
print(f"number of nodes with external_layer=False: {len(simulation.mesh.node_ids)}")

###############################################################################
# Extract stress/strain data
# --------------------------
# Extract stress or elastic strain data on the external layer.
# Averaging, and invariants computation can easily be done on the external layer since the
# connectivity of the external layer elements remains unchanged.

elemental_stress_ext = simulation.stress_principal_elemental(
    components=[1], external_layer=True
)
elemental_stress_ext.plot()

print(
    f"number of elements with external_layer=True: {len(elemental_stress_ext.index.mesh_index)}"
)
print(
    f"number of elements with external_layer=False: {len(simulation.mesh.element_ids)}"
)

###############################################################################
elastic_strain_eqv_ext = simulation.elastic_strain_eqv_von_mises_nodal(
    external_layer=True
)
elastic_strain_eqv_ext.plot()

###############################################################################
# Extract the external layer on a selection of elements
# -----------------------------------------------------

all_elements = simulation.mesh.element_ids
elements = []
for i in range(0, int(all_elements.size / 2)):
    elements.append(all_elements[i])
elemental_stress_ext = simulation.stress_principal_elemental(external_layer=elements)
elemental_stress_ext.plot()

###############################################################################
# Extract the external layer on a selection of elements for nodal results
# -----------------------------------------------------------------------

elastic_strain_eqv_ext = simulation.elastic_strain_eqv_von_mises_nodal(
    external_layer=elements
)
elastic_strain_eqv_ext.plot()

###############################################################################
# Extract the external layer on a selection of elements and scope results
# -----------------------------------------------------------------------

sub_elements = []
for i in range(0, int(len(elements) / 2)):
    sub_elements.append(elements[i])
elemental_stress_ext = simulation.stress_principal_elemental(
    external_layer=elements, element_ids=sub_elements
)
elemental_stress_ext.plot()
