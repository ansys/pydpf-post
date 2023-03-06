"""
.. _ref_invariants_example:

Extract stress/strain invariants - Static Simulation
=====================================================================
In this script a static simulation is used as an example to show how to
extract tensor's invariants.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.download_crankshaft()
simulation = post.load_simulation(example_path)

# for no autocompletion, this line is equivalent to:
simulation = post.StaticMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)

###############################################################################
# Extract elemental nodal stress and strain
# -----------------------------------------

stress = simulation.stress(all_sets=True)
print(stress)

strain = simulation.elastic_strain(all_sets=True)
print(strain)

###############################################################################
# Compute principal invariant averaged and unaveraged on stress and strain
# ------------------------------------------------------------------------

princ_stress_1 = simulation.stress_principal(components=[1])
print(princ_stress_1)

nodal_princ_stress_2 = simulation.stress_principal_nodal(components=[2])
print(nodal_princ_stress_2)
nodal_princ_stress_2.plot()

nodal_princ_strain_2 = simulation.elastic_strain_principal_nodal(components=[2])
print(nodal_princ_strain_2)
nodal_princ_strain_2.plot()


###############################################################################
# Compute Von Mises eqv averaged and unaveraged on stress and strain
# ------------------------------------------------------------------------

stress_eqv = simulation.stress_eqv_von_mises(set_ids=[1])
print(stress_eqv)

nodal_stress_eqv = simulation.stress_eqv_von_mises_nodal(set_ids=[1])
nodal_stress_eqv.plot()

nodal_strain_eqv = simulation.elastic_strain_eqv_von_mises_nodal(set_ids=[1])
nodal_strain_eqv.plot()
