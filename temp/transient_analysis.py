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

# Print the list of available named selections
named_selections = transient_simulation.named_selections
print(named_selections)

# Print time-steps
print(transient_simulation.time_freq_support)
print(transient_simulation.time_steps)

# Print the list of available results
print(transient_simulation.results)

# General plot of the simulation object with, by default:
# - the mesh, not deformed, at step 0
# - the geometry
# - the boundary conditions
transient_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)

# -----------------------------------------------------------------------------------------
# Apply a selection
selection = Selection()
selection.select_nodes(nodes=[1, 2, 3])
selection.select_time_freq_values(time_freq_values=[0.0, 0.1, 0.2])
transient_simulation.activate_selection(selection_object=selection)

# Deactivate a selection
transient_simulation.deactivate_selection()

# -----------------------------------------------------------------------------------------
# Extract results

# Extract displacements along X for nodes 1, 2 and 3 at t=0.05s
displacement_X = transient_simulation.displacement(
    component="X", nodes=[1, 2, 3], times=[0.05]
)
print(displacement_X)

# Extract nodal XY stresses for elements 1, 2 and 3 at time-step 1
stress_XY = transient_simulation.elemental_stress(
    component="XY", elements=[1, 2, 3], steps=[1]
)
print(stress_XY)

# Extract first principal nodal stress for a named (elemental or nodal) selection at all time-steps
stress_S1 = transient_simulation.nodal_stress(
    component="S1", named_selection=named_selections[0]
)
print(stress_S1)

# Extract equivalent elemental nodal strain for a selection at time-step 1
strain_eqv = transient_simulation.raw_strain(
    component="EQV", selection=selection, steps=[1]
)
print(strain_eqv)
