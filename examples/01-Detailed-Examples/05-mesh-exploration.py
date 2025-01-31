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

Explore the mesh
================
This example shows how to explore and manipulate the mesh object to query mesh data
such as connectivity tables, element IDs, element types and so on.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports.
# This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.common import elemental_properties

###############################################################################
# Load the result file
# --------------------
# Load the result file in a ``Simulation`` object that allows access to the results.
# The ``Simulation`` object must be instantiated with the path for the result file.
# For example, ``"C:/Users/user/my_result.rst"`` on Windows
# or ``"/home/user/my_result.rst"`` on Linux.

example_path = examples.download_harmonic_clamped_pipe()
simulation = post.HarmonicMechanicalSimulation(example_path)

###############################################################################
# Get the mesh
# ------------
# Retrieve the mesh and print it
mesh = simulation.mesh
print(mesh)

###############################################################################
# Plot the mesh
# -------------
# Plot the mesh to view the bare mesh of the model
mesh.plot()

###############################################################################
# Query basic information about the mesh
# --------------------------------------
# The ``Mesh`` object has several properties allowing access to different information such as:

###############################################################################
# the number of nodes
print(f"This mesh contains {mesh.num_nodes} nodes")

###############################################################################
# the list of node IDs
print(f"with IDs: {mesh.node_ids}")

###############################################################################
# the number of elements
print(f"This mesh contains {mesh.num_elements} elements")

###############################################################################
# the list of element IDs
print(f"with IDs {mesh.element_ids}")

###############################################################################
# the unit of the mesh
print(f"The mesh is in '{mesh.unit}'")

###############################################################################
# Named selections
# ----------------
# The available named selections are given as a dictionary
# with the names as keys and the actual ``NamedSelection`` objects as values.
# Printing the dictionary informs you on the available names.
named_selections = mesh.named_selections
print(named_selections)

###############################################################################
# To get a specific named selection, query it using its name as key
print(named_selections["_FIXEDSU"])

###############################################################################
# Elements
# --------
# Use ``mesh.elements`` to access the list of Element objects
print(mesh.elements)

###############################################################################
# You can then query a specific element by its ID
print(mesh.elements.by_id[1])

###############################################################################
# or by its index
element_0 = mesh.elements[0]
print(element_0)

###############################################################################
# Query information about a particular element
# --------------------------------------------
# You can request the IDs of the nodes attached to an ``Element`` object
print(element_0.node_ids)

###############################################################################
# or the list of ``Node`` objects
print(element_0.nodes)

###############################################################################
# To get the number of nodes attached, use
print(element_0.num_nodes)

###############################################################################
# Get the type of the element
print(element_0.type_info)
print(element_0.type)

###############################################################################
# Get the shape of the element
print(element_0.shape)

###############################################################################
# Element types and materials
# ---------------------------
# The ``Mesh`` object provides access to properties defined on all elements,
# such as their types or their associated materials.

###############################################################################
# Get the type of all elements
print(mesh.element_types)

###############################################################################
# Get the materials of all elements
print(mesh.materials)

###############################################################################
# Elemental connectivity
# ----------------------
# The elemental connectivity maps elements to connected nodes, either using IDs or indexes.

###############################################################################
# To access the indexes of the connected nodes using an element's index, use
element_to_node_connectivity = mesh.element_to_node_connectivity
print(element_to_node_connectivity[0])

###############################################################################
# To access the IDs of the connected nodes using an element's index, use
element_to_node_ids_connectivity = mesh.element_to_node_ids_connectivity
print(element_to_node_ids_connectivity[0])

###############################################################################
# Each connectivity object has a ``by_id`` property which changes the input from index to ID, thus:

###############################################################################
# To access the indexes of the connected nodes using an element's ID, use
element_to_node_connectivity_by_id = mesh.element_to_node_connectivity.by_id
print(element_to_node_connectivity_by_id[3487])

###############################################################################
# To access the IDs of the connected nodes using an element's ID, use
element_to_node_ids_connectivity_by_id = mesh.element_to_node_ids_connectivity.by_id
print(element_to_node_ids_connectivity_by_id[3487])

###############################################################################
# Nodes
# -----
# Get a node by its ID
node_1 = mesh.nodes.by_id[1]
print(node_1)

###############################################################################
# Get a node by its index
print(mesh.nodes[0])

###############################################################################
# Get the coordinates of all nodes
print(mesh.coordinates)

###############################################################################
# Query information about one particular node
# -------------------------------------------
# Get the coordinates of a node
print(node_1.coordinates)

###############################################################################
# Nodal connectivity
# ------------------
# The nodal connectivity maps nodes to connected elements, either using IDs or indexes.

###############################################################################
# To access the indexes of the connected elements using a node's index, use
node_to_element_connectivity = mesh.node_to_element_connectivity
print(node_to_element_connectivity[0])

###############################################################################
# To access the IDs of the connected elements using a node's index, use
node_to_element_ids_connectivity = mesh.node_to_element_ids_connectivity
print(node_to_element_ids_connectivity[0])

###############################################################################
# Each connectivity object has a ``by_id`` property which changes the input from index to ID, thus:

###############################################################################
# To access the indexes of the connected elements using a node's ID, use
node_to_element_connectivity_by_id = mesh.node_to_element_connectivity.by_id
print(node_to_element_connectivity_by_id[1])

###############################################################################
# To access the IDs of the connected elements using a node's ID, use
node_to_element_ids_connectivity_by_id = mesh.node_to_element_ids_connectivity.by_id
print(node_to_element_ids_connectivity_by_id[1])

###############################################################################
# Splitting into meshes
# ---------------------
# You can split the global mesh according to mesh properties to work on specific parts of the mesh
meshes = simulation.split_mesh_by_properties(
    properties=[elemental_properties.material, elemental_properties.element_shape]
)
###############################################################################
# The object obtained is a ``Meshes``
print(meshes)

###############################################################################
# Plotting a ``Meshes`` object plots a combination of all the ``Mesh`` objects within
meshes.plot(text="Mesh split")

###############################################################################
# Select a specific ``Mesh`` in the ``Meshes``, by index
meshes[0].plot(text="First mesh in the split mesh")

###############################################################################
# You can split the global mesh and select meshes based on specific property values
meshes_filtered = simulation.split_mesh_by_properties(
    properties={
        elemental_properties.material: [2, 3, 4],
        elemental_properties.element_shape: 1,
    }
)
meshes_filtered.plot(text="Mesh split and filtered")

###############################################################################
# or with a unique combination of property values
meshes[{"mat": 5, "elshape": 0}].plot(text="Mesh for mat=5 and elshape=0")
