"""Pseudocode for PyDPF-Post Transient Mechanical Simulations APIs."""

import ansys.dpf.post as dpf
from ansys.dpf.post import examples
from ansys.dpf.post.selection import Selection

# Provide Enums for available physics_types and analysis_types
# from ansys.dpf.post import physics_types, analysis_types


transient_simulation = dpf.load_simulation(examples.msup_transient)
# transient_simulation = dpf.load_simulation(examples.msup_transient,
#                                            physics_type=physics_types.mechanical,
#                                            analysis_type=analysis_types.transient)

# -----------------------------------------------------------------------------------------
# Explore the simulation metadata

# Print information about the simulation
print(transient_simulation)

# Print the mesh
print(transient_simulation.mesh)

# Print the list of constructed geometries
print(transient_simulation.geometries)

# Print the list of boundary conditions
print(transient_simulation.boundary_conditions)

# Print the list of loads
print(transient_simulation.loads)

# Print the list of available named selections
named_selections = transient_simulation.named_selections
print(named_selections)

# Print time-steps
print(transient_simulation.time_freq_support)
print(transient_simulation.time_steps)  # or steps or times

# Print the list of available results
print(transient_simulation.results)

# General plot of the simulation object with, by default:
# - the mesh, not deformed, at step 0
# - the geometry
# - the boundary conditions
transient_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)


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

# Mapping loads from one mesh to another
load = Load(data=data_object)
# A load would be derived from a DataObject, which has a mesh support, thus:
load.interpolate(mesh=another_mesh)


# -----------------------------------------------------------------------------------------
# Define and use a selection
selection = Selection()
selection.select_nodes(nodes=[1, 2, 3])
selection.select_elements(elements=[1, 2, 3])
selection.select_time_freq_indices(time_freq_indices=[1])  # Rework?
# selection.select_steps(steps=[1])
selection.select_named_selection(named_selection=named_selections[0])
selection.select_parts(part_ids=[1, 2])  # For LS-Dyna results.
# Could be available for all analyses types but would be undefined by default.


# -----------------------------------------------------------------------------------------
# Extract the mesh
mesh = transient_simulation.mesh

# Show the mesh with defaults for transient:
# - deformed if displacements are available
# - at last timestep
# - eroded elements not shown
mesh.plot(
    opacity=0.3, title="Transient mesh plot", text="defaults to deformed, last timestep"
)

# plot elements triads
mesh.plot(triads=True)


# -----------------------------------------------------------------------------------------
# Extract results

# Extract displacements along X for nodes 1, 2 and 3 at t=0.05s
displacement_X = transient_simulation.displacement(
    component="X", nodes=[1, 2, 3], times=[0.05]
)

# Extract norm of displacements for nodes 1, 2 and 3 at time-step 1
displacement_norm = transient_simulation.displacement(
    component="N", nodes=[1, 2, 3], steps=[1]
)

# Extract nodal XY stresses for elements 1, 2 and 3 at time-step 1
stress_XY = transient_simulation.elemental_stress(
    component="XY", elements=[1, 2, 3], steps=[1]
)

# Extract first principal nodal stress for a named (elemental or nodal) selection at all time-steps
stress_S1 = transient_simulation.nodal_stress(
    component="S1", named_selection=named_selections[0]
)

# Extract elemental Von Mises stress everywhere at all steps
stress_VM = transient_simulation.elemental_stress(component="VM")

# Extract equivalent elemental nodal elastic strain for a selection at step 1
elastic_strain_XY = transient_simulation.nodal_elastic_strain(
    component="XY", selection=selection, steps=[1]
)

# Extract first principal nodal strain for a selection at step 1
elastic_strain_E1 = transient_simulation.nodal_elastic_strain(
    component="E1", selection=selection, steps=[1]
)

# Extract nodal plastic strain for a selection at step 1
plastic_strain = transient_simulation.nodal_plastic_strain(
    selection=selection, steps=[1]
)

# Extract global internal energy at all times
global_internal_energy = transient_simulation.global_internal_energy()

# Extract hourglass energy for part 1 at all times
hourglass_energy_1 = transient_simulation.part_hourglass_energy(parts=[1])

# Extract hourglass energy for each part at all times
hourglass_energies = transient_simulation.part_hourglass_energy()

# -----------------------------------------------------------------------------------------
# Manipulate results with DataObject

# Plot a DataObject (3D contour plots)
# plot elemental hourglass energy at a given time on a specific deformed part
# argument "deformed" is True by default
elemental_hourglass_energy_1 = transient_simulation.elemental_hourglass_energy(
    parts=[1], times=[0.05]
)
elemental_hourglass_energy_1.plot(deformed=True)


# Create graphs using a DataObject (2D curve plots)
# plot an energy curve over time
global_internal_energy.plot()

# plot several energy curves over time
# Right now there is no logic to graph several curves using data from different DataObjects
# -> Need for combination of DataObjects (like for Dataframes)
# But then, what of the DataObject name? What about its support?
# -> Better to allow for curves from different DataObjects
# -> Use a GraphPlotter from Core,
# but create convenient helpers to quickly draw common energy graphs
