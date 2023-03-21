"""
.. _ref_mesh_exploration_example:

Explore the mesh
================
In this script a static simulation is used as an example to show how to
query mesh information such as connectivity, element IDs, element types and so on.
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

example_path = examples.download_crankshaft()
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

mesh_info = simulation.mesh_info  # TODO: expose MeshSelectionManager?
print(mesh_info)

exit()
###############################################################################
# Get the mesh
# ------------

# Use-case: big mesh, split by mat_id and thickness, get meshes for mat_id 1, limited evaluations

# print(simulation.mesh_info) -> Queries metadata
# ...
# Available properties/labels: "mat_id", "thickness", "eltype", "elshape"
# ...


# Mesh object holds on a MeshesContainer, specify splits and necessary evaluation defined at query
# Specific sub-mesh queries on Mesh by property

# Proposal 1: Mesh wraps a MeshesContainer - simulation.mesh getter
# meshes : Meshes = simulation.mesh(split_by="mat_id")
# meshes : Meshes = simulation.mesh(split_by=["mat_id", "thickness"])

# Get mesh split by thickness and mat_id but only for mat_id in 1, 2, 3
# mesh : Mesh = simulation.mesh(mat_id=[1, 2, 3], split_by=["mat_id", "thickness"])
# print(mesh)
# mesh : Mesh = simulation.mesh(mat_id=[1, 2, 3], thickness=)

# Pros:
# - Not lazy but very targeted evaluations -> same logic as result queries
# - simulation.mesh gets you a mesh/holder on meshes
# Cons:
# - simulation.mesh() by default would query everything as a single mesh?
# - not a property


# si


# simulation.mesh -> Mesh
# simulation.mesh.nodes_ids


# Proposal 2: Mesh wraps a MeshesContainer - no simulation.mesh, only query methods
# mesh_info : str = simulation.mesh_info

# Get the whole mesh
# mesh : Mesh = simulation.mesh (property)
# meshes : Meshes = mesh.split_by()

# Get the mesh split by material ID
# meshes : Meshes = simulation.split_mesh_by("mat_id", "thickness")  <--


# Get the mesh split by material ID but only for mat_id 1 and 2 and thickness 2
# meshes : Meshes = simulation.split_mesh_by(labels={"mat_id": [1, 2], "elshape": })  <--
# mesh = meshes.select(mat_id=1)
# simulation.create_elemental_named_selection() -> Selection
# simulation.create_named_selection(element_ids=[,2,,4], name="my_ns") -> Selection
# simulation.split_mesh_by({"named_selection"="my_ns")
# simulation.save_hfd5()
# mesh.name = "mat_id 1"


# print(meshes)
# """
# mesh1: {mat_id=1, thickness=2}
# mesh1: {mat_id=2, thickness=2}
# """


# print(meshes)
# """
# mesh0: {mat_id=1}
# mesh1: {mat_id=2}
# """
# Pros:
# - explicit
# Cons:


# Proposal 3: simulation.mesh is a "method" holder, not an object
# meshes : Mesh = simulation.mesh.get_mesh_split_by_property("mat_id", "thickness")
# print(meshes)
# """
# mesh0: {mat_id=1, thickness=1}
# mesh1: {mat_id=1, thickness=2}
# mesh2: {mat_id=2, thickness=1}
# mesh3: {mat_id=2, thickness=2}
# """

# Cons:
# simulation.mesh does not return anything


mesh = simulation.mesh
print(mesh)
# meshes : Meshes = simulation.meshes_with_elements_grouped_by("mat_id")

# meshes : Meshes = simulation.meshes(elements_grouped_by="mat_id")

# meshes : Meshes = simulation.meshes.elements_grouped_by("mat_id")

# meshes : Meshes = simulation.meshes_by_properties("mat_id", "thickness")
# print(meshes)
# """
# mesh0: {mat_id=1, thickness=1}
# mesh1: {mat_id=1, thickness=2}
# mesh2: {mat_id=2, thickness=1}
# mesh3: {mat_id=2, thickness=2}
# """
# mesh : Mesh = simulation.mesh_by_property(mat_id=1)
# mesh_info = simulation.mesh_info ->
# mesh = MeshesContainer.mesh_by_id(mat_id=1)
# for id, mesh for


###############################################################################
# Query basic information about the mesh (available)
# --------------------------------------

# Node IDs
n_ids = mesh.node_ids

# Element IDs
e_ids = mesh.element_ids

# Available named selection names
named_selections = mesh.available_named_selections

# Query basic information about the mesh (available in PyDPF-Core)
# --------------------------------------

# Number of nodes
n_nodes = mesh._core_object.nodes.n_nodes

# Number of elements
n_elements = mesh._core_object.elements.n_elements

# Number of named selections
n_ns = len(named_selections)

# Unit (get and set)
mesh_unit = mesh._core_object.unit
mesh._core_object.unit = "mm"

# Get the list of nodes/elements IDs of a given named selection
ns_scoping = mesh._core_object.named_selection(named_selection=named_selections[0])
ns_ids = ns_scoping.ids

# Named selection setter
mesh._core_object.set_named_selection_scoping(
    named_selection_name=named_selections[0],
    scoping=ns_scoping,
)

# Plot the mesh
plt = mesh._core_object.plot()

# #######################################
# Manipulating elements

# Get an element by ID
element_by_id = mesh._core_object.elements.element_by_id(1)
# Get an element by index
element_by_index = mesh._core_object.elements.element_by_index(0)
# Get the ID of an element
e_id = mesh._core_object.elements.element_by_index(0).id

# Adding elements to the mesh

# Get the element types
e_types = mesh._core_object.elements.element_types_field

# Get the materials
e_materials = mesh._core_object.elements.materials_field

# Get the elemental connectivity
connectivity = mesh._core_object.elements.connectivities_field

# ##############################################
# Query information about one particular element

# Get the nodes of an element
e_nodes = mesh._core_object.elements.element_by_id(1).nodes
# Get the node IDs of an element
e_node_ids = mesh._core_object.elements.element_by_id(1).node_ids
# Get the nodes of an element
e_n_node = mesh._core_object.elements.element_by_id(1).n_nodes

# Get the type of the element
e_type = mesh._core_object.elements.element_by_id(1).type
# Get the shape of the element
e_shape = mesh._core_object.elements.element_by_id(1).shape
# Get the connectivity of the element
e_connectivity = mesh._core_object.elements.element_by_id(1).connectivity


# ##################
# Manipulating nodes

# Get a node by ID
node_by_id = mesh._core_object.nodes.node_by_id(1)
# Get a node by index
node_by_index = mesh._core_object.nodes.node_by_index(0)

# Get the coordinates of all nodes
coordinates = mesh._core_object.nodes.coordinates_field

# ###########################################
# Query information about one particular node

# Coordinates
coor = mesh._core_object.nodes.node_by_id(1).coordinates
# Nodal connectivity
conn = mesh._core_object.nodes.node_by_id(1).nodal_connectivity
