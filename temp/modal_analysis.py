"""Pseudocode for PyDPF-Post Modal Mechanical Simulations APIs."""

import ansys.dpf.post as dpf
from ansys.dpf.post import examples
from ansys.dpf.post.common import AvailableSimulationTypes
from ansys.dpf.post.selection import Selection

# Load the simulation files
modal_simulation = dpf.load_simulation(
    examples.download_all_kinds_of_complexity_modal(),
    simulation_type=AvailableSimulationTypes.modal_mechanical,
)
# # Also possible:
# modal_simulation = dpf.load_simulation(examples.download_all_kinds_of_complexity_modal())
# modal_simulation = dpf.load_simulation(
#     examples.download_all_kinds_of_complexity_modal(), simulation_type="modal mechanical"
# )

# -----------------------------------------------------------------------------------------
# Explore the simulation metadata

# Print information about the simulation
print(modal_simulation)

# Print the mesh
print(modal_simulation.mesh)

# Print the list of constructed geometries
print(modal_simulation.constructed_geometries)

# Print the list of boundary conditions
print(modal_simulation.boundary_conditions)

# Print the list of available named selections
named_selections = modal_simulation.named_selections
print(named_selections)

# Print mode frequencies
print(modal_simulation.time_freq_support)
print(modal_simulation.modes)
print(modal_simulation.frequencies)

# Print the list of available results
print(modal_simulation.results)

# General plot of the simulation object with, by default:
# - the mesh, not deformed
# - the constructed geometry
# - the boundary conditions
modal_simulation.plot(mesh=True, constructed_geometries=True, boundary_conditions=True)

# -----------------------------------------------------------------------------------------
# Apply a selection
# Using the provided factories:
from ansys.dpf.post import tools

selection = tools.create_selection(
    node_ids=[1, 2, 3], element_ids=[1, 2, 3], load_steps=[1]
)
# or
selection = Selection(node_ids=[1, 2, 4], time_freq_indices=[0, 1])
selection = Selection()
selection.nodes(node_ids=[1, 2, 3])
selection.time_freq_indices(time_freq_indices=[0, 1])
modal_simulation.activate_selection(selection_object=selection)

# Deactivate a selection
modal_simulation.deactivate_selection()

# -----------------------------------------------------------------------------------------
# Extract results

# Extract displacements along X for nodes 1, 2 and 3 at f=0.05Hz
displacement_X = modal_simulation.displacement(
    components=["X"], nodes=[1, 2, 3], frequencies=[0.05]
)
print(displacement_X)

# Extract nodal XY stresses for elements 1, 2 and 3 at set 1
stress_XY = modal_simulation.elemental_stress(
    components=["XY"], elements=[1, 2, 3], set_ids=[1]
)
print(stress_XY)

# Extract first principal nodal stress for a named (elemental or nodal) selection at all frequencies
stress_S1 = modal_simulation.nodal_principal_stress(
    components=["1"], named_selection=named_selections[0]
)
print(stress_S1)

# Extract elemental equivalent strain for a selection at set 1
strain_eqv = modal_simulation.elemental_eqv_strain(selection=selection, set_ids=[1])
print(strain_eqv)

# -----------------------------------------------------------------------------------------
# Extract the mesh
mesh = modal_simulation.mesh

# Show the mesh with defaults for modal:
# - not deformed
# - at first mode/frequency
mesh.plot(opacity=0.3, title="Modal mesh plot", text="defaults to deformed, first mode")
