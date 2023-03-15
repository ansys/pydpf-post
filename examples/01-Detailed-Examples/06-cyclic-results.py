"""
.. _ref_cyclic_results_example:

Extract cyclic results
======================
In this script a modal analysis with cyclic symmetry is processed to show
how to expand the mesh and results.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf.core.plotter import DpfPlotter

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.download_multi_stage_cyclic_result()
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation.mesh._core_object)
# simulation.mesh._core_object.plot()
#############################################################################
# Extract expanded X displacement for the first mode
# --------------------------------------------------

ignore_cyclic = [True, False, False, False]
expand_sectors = [
    False,
    False,
]

for read_cyclic in range(4):
    displacement = simulation.displacement(components="X", read_cyclic=read_cyclic)
    fc = displacement._core_object
    f = fc[0]
    print("================")
    print(f"{read_cyclic=}")
    print(f"{displacement}")
    # print(f"{fc}")
    # print(f"{f}")
    # print(f"{f.meshed_region}")
    # print(f"{f.scoping.location=}")
    # print(f"{f.scoping.size=}")
    # print(f"{f.meshed_region.elements.n_elements=}")
    # print(f"{f.meshed_region.nodes.n_nodes=}")
    plt = DpfPlotter()
    try:
        plt.add_field(field=f)
    except:
        plt.add_field(field=f, meshed_region=simulation.mesh._core_object)
    plt.show_figure()
# print(displacement._core_object)
# print(displacement._core_object[0].meshed_region)  # Empty for read_cyclic=1
# displacement._core_object[0].meshed_region.plot()

# displacement.plot()
