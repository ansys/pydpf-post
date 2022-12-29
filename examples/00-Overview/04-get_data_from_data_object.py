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
print(f"Data length is {len(disp_last_time.ids)}.")

###############################################################################
# Get the displacement n_data
print(disp_last_time.n_data)

###############################################################################
# Get the displacement n_dim
print(disp_last_time.n_dim)

###############################################################################
# Get the displacement unit
disp_last_time.unit == "Pa"

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
print(disp_time25.name)
print(disp_time25.location)
print(disp_time25.unit)
print(disp_time25.data)
print(disp_time25.ids)
