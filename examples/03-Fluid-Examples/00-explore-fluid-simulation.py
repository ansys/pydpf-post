"""
.. _ref_explore_fluid_example:

Explore Fluid Simulation
========================
This example shows how to load a fluid simulation and explore the model and its available zones,
species, and phases.
"""
###############################################################################
# Use gRPC protocol on Linux
# --------------------------
# The CFF plugin is currently prone to errors when used InProcess on Linux,
# hence a gRPC server configuration is chosen when running on a Unix system.
import platform

import ansys.dpf.core as dpf
from ansys.dpf.core.dpf_operator import available_operator_names

from ansys.dpf import post
from ansys.dpf.post import examples

if "Linux" in platform.system():
    dpf.SERVER_CONFIGURATION = (
        dpf.server_factory.AvailableServerConfigs.LegacyGrpcServer
    )

# # # -> works for all results (no mach_number or mass_flow_rate)
# https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v231/en/cfx_tutr/i1361675.html
# fluid_example = examples.download_cfx_heating_coil()

# # -> works for all results (no mach_number or mass_flow_rate)
# https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v231/en/cfx_tutr/i1308373.html
fluid_example = examples.download_cfx_mixing_elbow()

# # -> does not work for mach_number, mass_flow_rate, surface_heat_rate, y_plus
# https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v231/en/flu_tg/flu_tg_elbow_mesh_wtc.html
# fluid_example = examples.download_fluent_mixing_elbow_steady_state()

# # # -> does not work for mach_number, mass_flow_rate, surface_heat_rate, y_plus
# https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v231/en/flu_tg/flu_tg_elbow_mesh_wtc.html
# fluid_example = examples.download_fluent_mixing_elbow_transient()

# # -> does not work for mass_flow_rate, surface_heat_rate
# https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v231/en/flu_tg/flu_tg_rotor_stator.html
# fluid_example = examples.download_fluent_axial_comp()
# opacity = 0.01

# # -> does not work for mach_number
# https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/corp/v231/en/flu_tg/flu_tg_magnus.html
# fluid_example = examples.download_fluent_multi_species()

server = dpf.SERVER
print(server.version)
available_operators = available_operator_names(server=server)
cff_operators = [operator for operator in available_operators if "cff::" in operator]
print("CFF operators =============================================")
print(cff_operators)
print("===========================================================")

# Load the fluid analysis result
cas = fluid_example["cas"]
dat = fluid_example["dat"]
simulation = post.FluidSimulation(cas=cas, dat=dat)
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

print(simulation.phases)
print(simulation.species)
print(simulation.zones)

for available_result in result_info.available_results:
    # if available_result.name in [
    #     "mach_number",
    #     "mass_flow_rate",
    #     "surface_heat_rate",
    #     "y_plus",
    # ]:
    #     continue
    print(available_result)
    result = getattr(simulation, available_result.name)()
    print(result)

    # result.plot()
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
