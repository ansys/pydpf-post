"""
.. _ref_explore_fluid_example:

Explore Fluid Simulation
========================
This example shows how to load a fluid simulation and explore the model and its available zones,
species, and phases.
"""
from ansys.dpf.core import examples
from ansys.dpf.core.dpf_operator import available_operator_names

from ansys.dpf import core as dpf
from ansys.dpf import post

fluid_example_files = examples.download_fluent_files()
ds = dpf.DataSources()
ds.set_result_file_path(
    fluid_example_files["cas"],
    key="cas",
)
ds.add_file_path(
    fluid_example_files["dat"],
    key="dat",
)

server = dpf.SERVER
print(server.version)
available_operators = available_operator_names(server=server)
cff_operators = [operator for operator in available_operators if "cff::" in operator]
print("CFF operators ==============================================")
print(cff_operators)
print("============================================================")

rip = dpf.operators.metadata.result_info_provider(data_sources=ds).eval()
print(rip)


# Load the fluid analysis result
simulation = post.FluidSimulation(ds)
print(simulation.result_info)

# # print(simulation.result_info.available_qualifier_labels)  # Requires new gate

# # density_op = model.operator("RHO")
# # density = density_op.outputs.fields_container()
density = simulation.density()
print(density)

density_on_nodes = simulation.density_on_nodes()
print(density_on_nodes)

density_on_faces = simulation.density_on_faces()
print(density_on_faces)

density_on_cells = simulation.density_on_cells()
print(density_on_cells)

# velocity = simulation.velocity()
# print(velocity)

# surface_heat_rate = simulation.surface_heat_rate()
# print(surface_heat_rate)
# # print(simulation.mesh_info)
# # print(simulation.zones)
# # print(simulation.species)
# # print(simulation.phases)
