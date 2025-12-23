# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
.. _ref_mesh_exploration_example:

Explore and manipulate the mesh
===============================
This example shows how to explore and manipulate the mesh to query mesh data
such as connectivity tables, element IDs, and element types.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.common import elemental_properties

###############################################################################
# Load result file
# ----------------
# Load the result file in a ``Simulation`` object that allows access to the results.
# The ``Simulation`` object must be instantiated with the path for the result file.
# For example, ``"C:/Users/user/my_result.rst"`` on Windows
# or ``"/home/user/my_result.rst"`` on Linux.

example_path = examples.download_harmonic_clamped_pipe()
simulation = post.HarmonicMechanicalSimulation(example_path)

###############################################################################
# Get mesh and print it
# ---------------------
mesh = simulation.mesh
print(mesh)

###############################################################################
# Plot mesh
# -------------
# Plot the mesh to view the bare mesh of the model.
mesh.plot()

###############################################################################
# Get basic information about mesh
# --------------------------------
# The ``Mesh`` object has several properties allowing access to different information.

###############################################################################
# Get the number of nodes.
print(f"This mesh contains {mesh.num_nodes} nodes")

###############################################################################
# Get the list of node IDs.
print(f"with IDs: {mesh.node_ids}")

###############################################################################
# Get the number of elements.
print(f"This mesh contains {mesh.num_elements} elements")

###############################################################################
# Get the list of element IDs.
print(f"with IDs {mesh.element_ids}")

###############################################################################
# Get the unit of the mesh.
print(f"The mesh is in '{mesh.unit}'")

###############################################################################
# Get named selections
# --------------------
# The available named selections are given as a dictionary
# with the names as keys and the actual ``NamedSelection`` objects as values.
#
# Print the dictionary to get the available names.
named_selections = mesh.named_selections
print(named_selections)

###############################################################################
# Get a specific named selection by using its name as the key.
print(named_selections["_FIXEDSU"])

###############################################################################
# Get elements
# --------
# Get a list of the elements.
print(mesh.elements)

###############################################################################
# Get a specific element by its ID.
print(mesh.elements.by_id[1])

###############################################################################
# Get a specific element by its index.
element_0 = mesh.elements[0]
print(element_0)

###############################################################################
# Get information about a particular element
# ------------------------------------------
# You can request the IDs of the nodes attached to an element.
print(element_0.node_ids)

###############################################################################
# Get the list of the element's nodes.
print(element_0.nodes)

###############################################################################
# Get the number of nodes attached to the element.
print(element_0.num_nodes)

###############################################################################
# Get the type of the element.
print(element_0.type_info)
print(element_0.type)

###############################################################################
# Get the shape of the element.
print(element_0.shape)

###############################################################################
# Get element types and materials
# -------------------------------
# The ``Mesh`` object provides access to properties defined on all elements,
# such as their types or associated materials.

###############################################################################
# Get the type of all elements.
print(mesh.element_types)

###############################################################################
# Get the materials of all elements.
print(mesh.materials)

###############################################################################
# Get elemental connectivity
# --------------------------
# The elemental connectivity maps elements to connected nodes using either IDs or indexes.

###############################################################################
# Access the indexes of the connected nodes using an element's index:
element_to_node_connectivity = mesh.element_to_node_connectivity
print(element_to_node_connectivity[0])

###############################################################################
# Access the IDs of the connected nodes using an element's index:
element_to_node_ids_connectivity = mesh.element_to_node_ids_connectivity
print(element_to_node_ids_connectivity[0])

###############################################################################
# Each connectivity object has a ``by_id`` property that changes the input from index to ID.

###############################################################################
# Access the indexes of the connected nodes using an element's ID.
element_to_node_connectivity_by_id = mesh.element_to_node_connectivity.by_id
print(element_to_node_connectivity_by_id[3487])

###############################################################################
# Access the IDs of the connected nodes using an element's ID:
element_to_node_ids_connectivity_by_id = mesh.element_to_node_ids_connectivity.by_id
print(element_to_node_ids_connectivity_by_id[3487])

###############################################################################
# Get a node or node information
# ------------------------------
# Get a node by its ID.
node_1 = mesh.nodes.by_id[1]
print(node_1)

###############################################################################
# Get a node by its index.
print(mesh.nodes[0])

###############################################################################
# Get the coordinates of all nodes.
print(mesh.coordinates)

###############################################################################
# Get the coordinates of a particular node.
print(node_1.coordinates)

###############################################################################
# Get nodal connectivity
# ----------------------
# The nodal connectivity maps nodes to connected elements, either using IDs or indexes.

###############################################################################
# Access the indexes of the connected elements using a node's index.
node_to_element_connectivity = mesh.node_to_element_connectivity
print(node_to_element_connectivity[0])

###############################################################################
# Access the IDs of the connected elements using a node's index.
node_to_element_ids_connectivity = mesh.node_to_element_ids_connectivity
print(node_to_element_ids_connectivity[0])

###############################################################################
# Each connectivity object has a ``by_id`` property that changes the input from index to ID.

###############################################################################
# Access the indexes of the connected elements using a node's ID.
node_to_element_connectivity_by_id = mesh.node_to_element_connectivity.by_id
print(node_to_element_connectivity_by_id[1])

###############################################################################
# Access the IDs of the connected elements using a node's ID.
node_to_element_ids_connectivity_by_id = mesh.node_to_element_ids_connectivity.by_id
print(node_to_element_ids_connectivity_by_id[1])

###############################################################################
# Split global mesh into mesh parts
# ---------------------------------
# You can split the global mesh according to mesh properties to work on specific parts of the mesh.
meshes = simulation.split_mesh_by_properties(
    properties=[elemental_properties.material, elemental_properties.element_shape]
)
###############################################################################
# A ``Meshes`` object obtained.
print(meshes)

###############################################################################
# Plot a ``Meshes`` object to plot a combination of all ``Mesh`` objects within the split mesh.
meshes.plot(text="Mesh split")

###############################################################################
# Select a specific ``Mesh``object in the split mesh by index.
meshes[0].plot(text="First mesh in the split mesh")

###############################################################################
# Split the global mesh and select meshes based on specific property values.
meshes_filtered = simulation.split_mesh_by_properties(
    properties={
        elemental_properties.material: [2, 3, 4],
        elemental_properties.element_shape: 1,
    }
)
meshes_filtered.plot(text="Mesh split and filtered")

###############################################################################
# Select a ``mesh`` object with a unique combination of property values.
meshes[{"mat": 5, "elshape": 0}].plot(text="Mesh for mat=5 and elshape=0")
