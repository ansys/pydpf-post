"""
.. _ref_fluid_example:

Fluid Simulation
================
This example shows how to load a fluid simulation and extract results like velocity and pressure.
It also shows extraction and selection of results by zone, phase or species.
"""

import ansys.dpf.core as dpf
from ansys.dpf.core.dpf_operator import available_operator_names

###############################################################################
# Perform required imports
# ------------------------
from ansys.dpf import post
from ansys.dpf.post import examples

# Start a server and check its version
server = dpf.start_local_server()
print(f"{server.version=}")
# # Check CFF operators are now available
# available_operators = available_operator_names(server=server)
# cff_operators = [operator for operator in available_operators if "cff::" in operator]
# print("CFF operators ==============================================")
# print(cff_operators)
# print("============================================================")

# Load the fluid analysis result
simulation = post.FluidSimulation(
    r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj"  # noqa
)
print(simulation)

mesh = simulation.mesh
print(mesh)
# mesh.plot()

available_operators = available_operator_names(server=server)
cff_operators = [operator for operator in available_operators if "cff::" in operator]
print("CFF operators ==============================================")
print(cff_operators)
print("============================================================")

# Explore the structure of the results
# print(simulation.zones)
# print(simulation.species)
# print(simulation.phases)

# Extract a result
density = simulation.density(zone_ids=[-1], phases=[1])
print(density)
