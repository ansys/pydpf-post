"""Pseudocode for PyDPF-Post Modal Mechanical Simulations APIs."""

import ansys.dpf.post as dpf
from ansys.dpf.post import examples
from ansys.dpf.post.selection import Selection

# Provide Enums for available physics_types and analysis_types
# from ansys.dpf.post import physics_types, analysis_types


modal_simulation = dpf.load_simulation(
    examples.download_all_kinds_of_complexity_modal()
)
# modal_simulation = dpf.load_simulation(examples.download_all_kinds_of_complexity_modal(),
#                                        physics_type=physics_types.mechanical,
#                                        analysis_type=analysis_types.modal)

# -----------------------------------------------------------------------------------------
# Explore the simulation metadata

# Print information about the simulation
print(modal_simulation)

# Print the mesh
print(modal_simulation.mesh)

# Print the list of constructed geometries
print(modal_simulation.geometries)

# Print the list of boundary conditions
print(modal_simulation.boundary_conditions)

# Print the list of available named selections
named_selections = modal_simulation.named_selections
print(named_selections)

# Print mode frequencies
print(modal_simulation.time_freq_support)

# Print the list of available results
print(modal_simulation.results)

# General plot of the simulation object with, by default:
# - the mesh, not deformed
# - the geometry
# - the boundary conditions
modal_simulation.plot(mesh=True, geometry=True, boundary_conditions=True)

# -----------------------------------------------------------------------------------------
# Apply a selection
selection = Selection()
selection.select_nodes(nodes=[1, 2, 3])
selection.select_time_freq_indices(time_freq_indices=[0, 1])
modal_simulation.activate_selection(selection_object=selection)

# Deactivate a selection
modal_simulation.deactivate_selection()

# -----------------------------------------------------------------------------------------
# Extract results

# Extract displacements along X for nodes 1, 2 and 3 at f=0.05Hz
displacement_X = modal_simulation.displacement(
    component="X", nodes=[1, 2, 3], frequencies=[0.05]
)
print(displacement_X)

# Extract nodal XY stresses for elements 1, 2 and 3 at frequency 1
stress_XY = modal_simulation.elemental_stress(
    component="XY", elements=[1, 2, 3], steps=[1]
)
print(stress_XY)

# Extract first principal nodal stress for a named (elemental or nodal) selection at all frequencies
stress_S1 = modal_simulation.nodal_stress(
    component="S1", named_selection=named_selections[0]
)
print(stress_S1)

# Extract equivalent elemental nodal strain for a selection at the second frequency
strain_eqv = modal_simulation.raw_strain(
    component="EQV", selection=selection, steps=[1]
)
print(strain_eqv)
