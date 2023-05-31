"""
.. _ref_fluid_example:

Fluid Simulation
================
This example shows how to load a fluid simulation and extract results like velocity and pressure.
It also shows extraction and selection of results by zone, phase or species.
"""
###############################################################################
# Perform required imports
# ------------------------
from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.post import examples

# Load the fluid analysis result
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

simulation = post.FluidSimulation(ds)
print(simulation)

# mesh = simulation.mesh
# print(mesh)
# mesh.plot()

# Explore the structure of the results
# print(simulation.zones)
# print(simulation.species)
# print(simulation.phases)

# Extract a result
density = simulation.density(zone_ids=[-1], phases=[1])
print(density)
