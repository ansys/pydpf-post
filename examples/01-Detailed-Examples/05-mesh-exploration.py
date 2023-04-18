"""
.. _ref_mesh_exploration_example:

Explore the mesh
================
In this script a static simulation is used as an example to show how to
query mesh information such as connectivity, element IDs, element types and so on.
"""
from __future__ import annotations

from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.common import elemental_properties

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.


###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.download_all_kinds_of_complexity()
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)
stress_df = simulation.stress(components=["XX"])
print(stress_df)

###############################################################################
# Get the mesh
# ------------

# Retrieve the actual mesh
mesh = simulation.mesh

###############################################################################
# Query basic information about the mesh (available in PyDPF-Core)
# --------------------------------------

# Node IDs
n_ids = mesh.node_ids

# Element IDs
e_ids = mesh.element_ids

# Available named selection names
named_selections = mesh.available_named_selections

# Number of nodes
n_nodes = mesh.num_nodes

# Number of elements
n_elements = mesh.num_elements

# Number of named selections
n_ns = len(named_selections)

# Unit (get and set)
mesh_unit = mesh.unit
mesh.unit = "mm"

print(n_ids)
print(e_ids)
print(named_selections)
print(n_nodes)
print(n_elements)

###############################################################################
# Get Named Selections
# --------------------
ns_list     = mesh.available_named_selections
first_key = ns_list[0]
named_selection = mesh.named_selections[first_key]

for k in mesh.named_selections.keys():
    print(k)
for v in mesh.named_selections.values():
    print(v)
for k, v in mesh.named_selections.items():
    print(f"{k} = {v}")

###############################################################################
# Get elements
# ------------

# Get an element by ID
el_by_id = mesh.elements.by_id[1]

# Get an element by index
index = el_by_id.index
print(mesh.elements[index])

###############################################################################
# Element Types and Materials
# ---------------------------

# Get the element types
el_types = mesh.element_types
print(el_types)

# Get the materials
e_materials = mesh.materials
print(e_materials)

###############################################################################
# Query information about one particular element
# ----------------------------------------------

# Get the nodes of an element
mesh.elements[1].nodes

# Get the node IDs of an element
mesh.elements[1].node_ids

# Get the nodes of an element
mesh.elements[1].n_nodes

# Get the type of the element
mesh.elements[1].type_info
mesh.elements[1].type_id

# Get the shape of the element
mesh.elements[1].shape

###############################################################################
# Get the elemental connectivity
# ------------------------------

# get node indices from element index
conn1 = mesh.conn_elem_to_node

# get node IDs from element index
conn2 = mesh.conn_elem_to_node_id

el_idx_5 = mesh.elements[5]
# get node IDS from element ID
conn2.by_id[el_idx_5.id]

###############################################################################
# Get nodes
# ---------

# Get a node by ID
node_by_id = mesh.nodes.by_id[1]

# Get a node by index
node_by_index = mesh.nodes[0]

# Get the coordinates of all nodes
print(mesh.coordinates)

###############################################################################
# Query information about one particular node
# -------------------------------------------

# Coordinates
mesh.nodes[1].coordinates

# Get Nodal connectivity
conn3 = mesh.conn_node_to_elem
conn4 = mesh.conn_node_to_elem_id

print(mesh.nodes[0])

# elements indices to elem_id
print(conn3[0])

# elements IDs from node index
print(conn4[0])

###############################################################################
# Splitting into meshes
# --------

# Plot the mesh
mesh.plot()

# Split the global mesh according to mesh properties
meshes = simulation.split_mesh_by_properties(
    properties=[elemental_properties.material, elemental_properties.element_shape]
)
meshes.plot()

# Split the global mesh and select meshes for specific property values
print(meshes)
meshes = simulation.split_mesh_by_properties(
    properties={
        elemental_properties.material: 1,
        elemental_properties.element_shape: [0, 1],
    }
)

# Mesh<index=0, num_nodes=100, num_elem=1000>
meshes.plot()

# Select a specific Mesh in the Meshes, by index
meshes[1].plot()
# or by property values
meshes[{"mat": 1, "elshape": 0}].plot()


