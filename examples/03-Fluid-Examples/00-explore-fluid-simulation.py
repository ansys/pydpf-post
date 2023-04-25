"""
.. _ref_explore_fluid_example:

Explore Fluid Simulation
========================
This example shows how to load a fluid simulation and explore the model and its available zones,
species, and phases.
"""
from ansys.dpf import post

# Load the fluid analysis result
simulation = post.FluidSimulation(
    r"D:\ANSYSDev\Sandbox\plugins\Ans.Dpf.CFF\source\Ans.Dpf.CFFTest\test_models\FLPRJ\axial_comp\axial_comp_reduced.flprj"  # noqa
)
print(simulation)

print(simulation.mesh)
# print(simulation.zones)
# print(simulation.species)
# print(simulation.phases)
