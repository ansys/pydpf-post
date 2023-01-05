"""Pseudocode for PyDPF-Post Static Mechanical Simulations APIs."""

import ansys.dpf.post as dpf
from ansys.dpf.post import examples
from ansys.dpf.post.selection import Selection

# Provide Enums for available physics_types and analysis_types
# from ansys.dpf.post import physics_types, analysis_types

# Load the simulation files
static_simulation = dpf.load_simulation(examples.simple_bar)
# static_simulation = dpf.load_simulation(examples.simple_bar,
#                                         physics_type=physics_types.mechanical,
#                                         analysis_type=analysis_types.static)

# -----------------------------------------------------------------------------------------
# Explore the simulation metadata

# Print information about the simulation
print(static_simulation)

# Print the mesh
print(static_simulation.mesh)

# Print the list of constructed geometries
print(static_simulation.geometries)

# Print the list of boundary conditions
print(static_simulation.boundary_conditions)

# Print the list of loads
print(static_simulation.loads)

# Print the list of available named selections
named_selections = static_simulation.named_selections
print(named_selections)

# Print the list of steps (the time_freq_support)
print(static_simulation.time_freq_support)
print(static_simulation.steps)

# Print the list of available results
print(static_simulation.results)

# General plot of the simulation object with, by default:
# - the mesh, not deformed, at step 0
# - the geometry
# - the boundary conditions
static_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)


# -----------------------------------------------------------------------------------------
# Create and use geometry
# Creation of geometry is for now in PyDPF-Core.
# Do we expose those in Post in some way? Simply by importing them in a specific Post module?


# -----------------------------------------------------------------------------------------
# Define and use boundary conditions
# TODO: TBD


# -----------------------------------------------------------------------------------------
# Define and use loads
# TODO: TBD

# Loads are not read from files yet
# A load would then be derived from a DataObject, which has a mesh support, thus:
load = Load(data=data_object)
# Mapping loads from one mesh to another (which is the same as DataObject.interpolate?)
load.interpolate(mesh=another_mesh)
# This means a load has a support, either nodal, elemental or face.
# This means location.face must be available.
# It can also vary with step/time/frequency/rpm or any dimension present in the DataObject.

# The metadata and functionalities unique to a Load would be:
# - its naming
# - its string representation
# - its graphic representation

# -----------------------------------------------------------------------------------------
# Define and use a selection
selection = Selection()
selection.select_nodes(nodes=[1, 2, 3])
selection.select_elements(elements=[1, 2, 3])
selection.select_time_freq_indices(time_freq_indices=[1])  # Rework?
# selection.select_steps(steps=[1])
selection.select_geometry(geometry=Line())
selection.select_named_selection(named_selection=named_selections[0])

# Intersect two spatial selections
selection_1 = Selection()
selection_1.select_nodes(nodes=[1, 2, 3])
selection_2 = Selection()
selection_2.select_nodes(nodes=[3, 4, 5])

space_sel_1 = selection_1.spatial_selection
print(space_sel_1)
# Spatial selection with 3 node(s):
# [1, 2, 3]
space_sel_2 = selection_2.spatial_selection
print(space_sel_2)
# Spatial selection with 3 node(s):
# [3, 4, 5]

space_sel_1.intersect(space_sel_2)
print(space_sel_1)
# Spatial selection with 1 node(s):
# [3]

# Union of two spatial selections
selection_3 = Selection()
selection_3.select_nodes(nodes=[1, 2, 3])
selection_4 = Selection()
selection_4.select_nodes(nodes=[3, 4, 5])
space_sel_3 = selection_3.spatial_selection
space_sel_4 = selection_4.spatial_selection

space_sel_3.union(space_sel_4)
print(space_sel_3)
# Spatial selection with 5 node(s):
# [1, 2, 3, 4, 5]

# Do we then need to explicitly update the spatial_selection of selection_1?

# Select based on constructed geometry
# TODO

# Set a selection as default for a simulation for further result extractions
static_simulation.activate_selection(selection_object=selection)

# Remove a selection as default for a simulation for further result extractions
static_simulation.deactivate_selection()

# Apply a SpatialSelection to a Simulation to get resulting list of entities IDs (nodes or elements)
list_of_IDs = selection.spatial_selection.apply_to(static_simulation)

# -----------------------------------------------------------------------------------------
# Extract the mesh
mesh = static_simulation.mesh
print(mesh)

# Show the mesh with defaults for static:
# - not deformed
# - at last step
mesh.plot(
    opacity=0.3, title="Static mesh plot", text="defaults to not deformed, last step"
)

# plot elements triads
mesh.plot(triads=True)


# -----------------------------------------------------------------------------------------
# Extract results

# Extract displacements along X for nodes 1, 2 and 3 at step 1
displacement_X = static_simulation.displacement(
    component="X", nodes=[1, 2, 3], steps=[1]
)

# Extract norm of displacements for nodes 1, 2 and 3 at step 1
displacement_norm = static_simulation.displacement(
    component="N", nodes=[1, 2, 3], steps=[1]
)

# Extract nodal XY stresses for elements 1, 2 and 3 at step 1
stress_XY = static_simulation.elemental_stress(
    component="XY", elements=[1, 2, 3], steps=[1]
)

# Extract first principal nodal stress for a named (elemental or nodal) selection at all steps
stress_S1 = static_simulation.nodal_stress(
    component="S1", named_selection=named_selections[0]
)

# Extract elemental Von Mises stress everywhere at all steps
stress_VM = static_simulation.elemental_stress(component="VM")

# Extract equivalent elemental nodal elastic strain for a selection at step 1
elastic_strain_XY = static_simulation.nodal_elastic_strain(
    component="XY", selection=selection, steps=[1]
)

# Extract first principal nodal strain for a selection at step 1
elastic_strain_E1 = static_simulation.nodal_elastic_strain(
    component="E1", selection=selection, steps=[1]
)

# Extract nodal plastic strain for a selection at step 1
plastic_strain = static_simulation.nodal_plastic_strain(selection=selection, steps=[1])

# -----------------------------------------------------------------------------------------
# Manipulate results with DataObject

# Print a DataObject
print(stress_S1)
# Nodal stress S1
# For:
# - steps: 1, 2, 3, ...
# - named selection: {named_selection[0]}

# Get the name of the DataObject
print(stress_S1.name)
# Nodal stress S1

# Get the mesh support of the DataObject
print(stress_S1.mesh)

# Get the shape of the DataObject
print(stress_S1.shape)
# list of lengths of each dimension present: [nb_steps, nb_entities,

# Get specific data from a DataObject (as a new DataObject)
stress_S1_step_1 = stress_S1[1]  # Equivalent to stress_S1[0, ...] or stress_S1[0, :, :]
# or
# stress_S1_step_1 = stress_S1.step(1)
print(stress_S1_step_1)
# Nodal stress S1
# For:
# - steps: 1
# - named selection: {named_selection[0]}


# Get the minimum, maximum, mean of the DataObject (using DPF operators for efficiency)
# Across all dimensions present in the specific DataObject, unless specified as argument
# Take inspiration from the "axis" argument for ``numpy.ndarray.amax``.
# Do we want a ``nanmax``? Which propagates NaN values? (Whether we have actual NaNs or not)
# These would return DataObjects or a scalar.
global_max_stress_S1 = stress_S1.max()

max_stress_S1_step_1 = stress_S1_step_1.max()
# or
# max_stress_S1_step_1 = stress_S1.max(steps=[1])

# Arguments would include lists of steps, nodes, elements, components... a named selection
max_stress_S1_node_1 = stress_S1.max(nodes=[1])


# Get the mesh support for a specific DataObject (uses the scoping associated to the DataObject)
mesh_1 = max_stress_S1_node_1.mesh
print(mesh_1)
# Mesh containing:
# - 1 nodes: [1]
# - 0 elements: []


# Plot a DataObject (3D contour plots)
# TODO
# plotting and graphing functions can return the Plotter instance
# plot a multi-step DataObject - takes an optional "step" argument, defaults to last
pl = stress_S1.plot(step=3, return_plotter=True)

# plot a one-step DataObject - the "step" argument is not read
stress_S1_step_1.plot()

# plot principal stresses
# TODO: Not really sure how this is different from plotting results


# Create graphs using a DataObject (2D curve plots)
# TODO
# the graph method would take optional "x", "y", "names" arguments.
# graph a multi-step DataObject
max_stress_S1_node_1.graph(x="steps", names="ID")
# or as proposed in the GitHub discussion #213
max_stress_S1_node_1.plot(graph=True, individual_figures=False)


# Animate a DataObject
# animate a contour in time
stress_S1.animate(axis="steps")


# Animate graphs?
# animate a curve in time
max_stress_S1_node_1.animate_graph(axis="steps", x="ID")


# Interpolate data from one mesh to another
# Starts from the mesh support of the DataObject
# Performs interpolation for each step and for each component present by default (for each field)
stress_S1_on_another_mesh = stress_S1.interpolate(mesh=another_mesh)


# Export to a numpy ndarray
displacement_X_arr = displacement_X.as_array()

# Export to a pandas Dataframe
displacement_X_df = displacement_X.as_data_frame()
