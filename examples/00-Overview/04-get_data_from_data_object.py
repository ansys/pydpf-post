"""
.. _ref_get_data_from_data_object:

Get data from DataObject
========================
This example show how to interact with the ``DataObject`` loaded from the simulation.
"""

###############################################################################
# Imports and loading simulation
# ------------------------------
from ansys.dpf.post import examples, load_simulation

simulation = load_simulation(examples.download_transient_result())

###############################################################################
# Get displacements from simulation
# ---------------------------------
# Get the displacements for all time steps
dispObject = simulation.displacement()
print(dispObject)

###############################################################################
# Get the one entry of the displacements. Select the last time step
disp_time0 = dispObject[-1]
print(disp_time0)

###############################################################################
# Get the displacement name
disp_last_time = dispObject[-1]
print(disp_last_time.name)

###############################################################################
# Get the displacement location
print(disp_last_time.location)

###############################################################################
# Get the displacement data
print(disp_last_time.data)
print(f"Data length is {len(disp_last_time.data)}.")

###############################################################################
# Get the displacement IDs
print(disp_last_time.ids)
print(f"IDs length is {len(disp_last_time.ids)}.")

###############################################################################
# Get the displacement n_data
print(disp_last_time.n_data)

###############################################################################
# Get the displacement n_dim
print(disp_last_time.n_dim)

###############################################################################
# Get the displacement unit
print(disp_last_time.unit)

###############################################################################
# Plot displacement field
disp_last_time.plot()

###############################################################################
# Get displacements for a single time step and a few nodes
# --------------------------------------------------------
dispObject = simulation.displacement(nodes=[3, 4, 6], steps=[1, 25])
print(dispObject)

###############################################################################
# Inspect one of the two entries of ``dispObject``
disp_time25 = dispObject[1]
print(f"Name = {disp_time25.name}")
print(f"Location = {disp_time25.location}")
print(f"Unit = {disp_time25.unit}")
print(f"Data = {disp_time25.data}")
print(f"IDs = {disp_time25.ids}")
print(f"N_data = {disp_time25.n_data}")
print(f"N_dim = {disp_time25.n_dim}")

###############################################################################
# Get elemental stresses for one component and a few elements
# --------------------------------------------------------
stressObject = simulation.elemental_stress(elements=[14, 22, 75], component="XY")
print(stressObject)

###############################################################################
# Inspect one entry of ``stressObject``
stress_xy = stressObject[0]
print(f"Name = {stress_xy.name}")
print(f"Location = {stress_xy.location}")
print(f"Unit = {stress_xy.unit}")
print(f"Data = {stress_xy.data}")
print(f"IDs = {stress_xy.ids}")
print(f"N_data = {stress_xy.n_data}")
print(f"N_dim = {stress_xy.n_dim}")
