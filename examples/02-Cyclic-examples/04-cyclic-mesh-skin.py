"""
.. _ref_mesh_skin_example:

Extract Data and Plot on Mesh Skin
==================================
In this example, post-processing on a mesh skin for a static analysis is displayed.
The skin mesh is rebuilt with new surface elements connecting the nodes on the external skin
of the solid mesh. These surface elements types are chosen with respect to the solid elements
facets having all their nodes on the skin.
This feature is available for all types of MechanicalSimulation and allows you to reduce the size
of the mesh and of the extracted data to improve processing performances. Larger stress and strains
being usually located on the skin of a model, computing results on the external layer gives
equivalent maximum values in most cases. Post-processing of elemental or elemental nodal results
requires element solid to skin mapping in order to get from a solid element result to a facet 
result.
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
# Extract displacement data over the entire mesh and on the external layer.

displacement_skin = simulation.displacement(skin=True)
displacement = simulation.displacement()  # default is skin=False
displacement_skin.plot()
displacement.plot()

print(
    f"number of nodes with skin=True: {len(displacement_skin.index.mesh_index)}"
)
print(
    f"number of nodes with skin=False: {len(displacement.index.mesh_index)}"
)

###############################################################################
# Extract stress/strain data
# --------------------------
# Extract stress or elastic strain data over the entire mesh and on the external layer.
# Averaging, and invariants computation can easily be done on the external layer since the
# connectivity of the kept elements remains unchanged.

elemental_stress_skin = simulation.stress_principal_elemental(components=[1], skin=True)
elemental_stress = simulation.stress_principal_elemental()
elemental_stress_skin.plot()
elemental_stress.plot()

print(
    f"number of elements with skin=True: {len(elemental_stress_skin.index.mesh_index)}"
)
print(
    f"number of elements with skin=False: {len(elemental_stress.index.mesh_index)}"
)


elastic_strain_eqv_skin = simulation.elastic_strain_eqv_von_mises_nodal(skin=True)
elastic_strain_eqv = simulation.elastic_strain_eqv_von_mises_nodal()
elastic_strain_eqv_skin.plot()
elastic_strain_eqv.plot()
