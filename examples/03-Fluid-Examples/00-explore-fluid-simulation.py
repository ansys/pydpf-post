"""
.. _ref_explore_fluid_example:

Explore Fluid Simulation
========================
This example shows how to load a fluid simulation and explore the model and its available zones,
species, and phases.
"""
from ansys.dpf.core.dpf_operator import available_operator_names

from ansys.dpf import post

# Load the fluid analysis result
simulation = post.FluidSimulation(
    r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj"  # noqa
)
server = simulation._model._server
available_operators = available_operator_names(server=server)
cff_operators = [operator for operator in available_operators if "cff::" in operator]
print("CFF operators ==============================================")
print(cff_operators)
print("============================================================")

print(simulation.result_info)
# print(simulation.result_info.available_qualifier_labels)
print(simulation.density())

# print(simulation.mesh)
# print(simulation.zones)
# print(simulation.species)
# print(simulation.phases)
