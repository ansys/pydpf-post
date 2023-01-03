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

# Print the list of available named selections
named_selections = static_simulation.named_selections
print(named_selections)

# Print the list of steps (the time_freq_support)
print(static_simulation.time_freq_support)

# Print the list of available results
print(static_simulation.results)

# General plot of the simulation object with, by default:
# - the mesh, not deformed, at step 0
# - the geometry
# - the boundary conditions
static_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)

# -----------------------------------------------------------------------------------------
# Apply a selection
selection = Selection()
selection.select_nodes(nodes=[1, 2, 3])
selection.select_time_freq_indices(time_freq_indices=[1])
static_simulation.activate_selection(selection_object=selection)

# Deactivate a selection
static_simulation.deactivate_selection()

# -----------------------------------------------------------------------------------------
# Extract results

# Extract displacements along X for nodes 1, 2 and 3 at step 1
displacement_X = static_simulation.displacement(
    component="X", nodes=[1, 2, 3], steps=[1]
)
print(displacement_X)

# Extract nodal XY stresses for elements 1, 2 and 3 at step 1
stress_XY = static_simulation.elemental_stress(
    component="XY", elements=[1, 2, 3], steps=[1]
)
print(stress_XY)

# Extract first principal nodal stress for a named (elemental or nodal) selection at all steps
stress_S1 = static_simulation.nodal_stress(
    component="S1", named_selection=named_selections[0]
)
print(stress_S1)

# Extract equivalent elemental nodal strain for a selection at step 1
strain_eqv = static_simulation.raw_strain(
    component="EQV", selection=selection, steps=[1]
)
print(strain_eqv)
