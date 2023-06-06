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

# fluid_example = examples.download_fluent_multi_species()  # -> mach_number
# ds = dpf.DataSources()
# ds.set_result_file_path(fluid_example["cas"], "cas")
# ds.add_file_path(fluid_example["dat"], "dat")

# # -> works for all results (no mach_number or mass_flow_rate)
# fluid_example = examples.download_cfx_heating_coil()
# ds = dpf.DataSources()
# ds.set_result_file_path(fluid_example["cas"], "cas")
# ds.add_file_path(fluid_example["dat"], "dat")

# # -> works for all results (no mach_number or mass_flow_rate)
# fluid_example = examples.download_cfx_mixing_elbow()
# ds = dpf.DataSources()
# ds.set_result_file_path(fluid_example["cas"], "cas")
# ds.add_file_path(fluid_example["dat"], "dat")

# fluid_example = examples.download_fluent_mixing_elbow_steady_state()  # -> mach_number
# ds = dpf.DataSources()
# ds.set_result_file_path(fluid_example["cas"][0], "cas")
# ds.add_file_path(fluid_example["dat"][0], "dat")

fluid_example = examples.download_fluent_mixing_elbow_transient()  # -> mach_number
ds = dpf.DataSources()
ds.set_result_file_path(fluid_example["cas"][0], "cas")
ds.add_file_path(fluid_example["dat"][0], "dat")

# fluid_example = examples.download_fluent_axial_comp()  # -> mass_flow_rate fails
# print(fluid_example)
# ds = dpf.DataSources()
# ds.set_result_file_path(fluid_example["cas"][0], "cas")
# ds.add_file_path(fluid_example["dat"][0], "dat")

# ds.set_result_file_path(
#     fluid_example_files["cas"],
#     key="cas",
# )
# ds.add_file_path(
#     fluid_example_files["dat"],
#     key="dat",
# )

server = dpf.SERVER
print(server.version)
available_operators = available_operator_names(server=server)
cff_operators = [operator for operator in available_operators if "cff::" in operator]
print("CFF operators ==============================================")
print(cff_operators)
print("============================================================")

# Load the fluid analysis result
simulation = post.FluidSimulation(ds)
print(simulation)

# print(simulation.results)
# print(simulation.result_info.available_results[3])
# mach_op = simulation._model.operator("MACH")
# # mach_op.connect(25, 5)
# mach_number = mach_op.eval()
# print(mach_number)
#
# exit()

result_info = simulation.result_info
# print(result_info)

# print(simulation.phases)
# print(simulation.species)

for available_result in result_info.available_results:
    print(available_result)
    print(getattr(simulation, available_result.name)(zone_ids=[1]))

# getattr(simulation, available_result.name)().plot(opacity=0.5)
# print(dir(result_info))
# print(result_info.available_qualifier_labels)
# print(result_info.qualifier_label_support)
# simulation.mesh._meshed_region.plot(opacity=0.3)

# # print(simulation.result_info.available_qualifier_labels)  # Requires new gate

# model = simulation._model

# density_op = model.operator("RHO")
# density = density_op.outputs.fields_container()
# print(density)
# density[0].plot()

# for result_name in simulation.result_info.available_results:
#     result = getattr(simulation, result_name.name)()
#     print(result)
#     result.plot(opacity=0.5)

# mean_velocity = simulation.mean_velocity()
# print(mean_velocity)
# mean_velocity.plot(opacity=0.5)
#
# density_on_nodes = simulation.density_on_nodes()
# print(density_on_nodes)
#
# density_on_faces = simulation.density_on_faces()
# print(density_on_faces)
#
# density_on_cells = simulation.density_on_cells()
# print(density_on_cells)

# velocity = simulation.velocity()
# print(velocity)

# surface_heat_rate = simulation.surface_heat_rate()
# print(surface_heat_rate)
# # print(simulation.mesh_info)
# # print(simulation.zones)
# # print(simulation.species)
# # print(simulation.phases)
