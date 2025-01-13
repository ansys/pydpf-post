"""
.. _ref_basics:

Get and use a ``Result`` object
===============================
This example shows how to use the legacy PyDPF-Post API to get and use a ``Result`` object.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Solution`` object
# -----------------------
# Get the ``Solution`` object and instantiate with the path for the result
# file.

example_path = examples.download_all_kinds_of_complexity()
solution = post.load_solution(example_path)

###############################################################################
# Get ``Result`` objects
# ----------------------

###############################################################################
# Get displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~

displacement_result = solution.displacement()
displacement = displacement_result.vector

###############################################################################
# Get information on displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

displacement.num_fields
disp_data = displacement.get_data_at_field(0)
len(disp_data)

disp_data[1]

displacement.max_data
displacement.get_max_data_at_field(0)

displacement.min_data

###############################################################################
# Get tensor stress result
# ~~~~~~~~~~~~~~~~~~~~~~~~
# You can get the nodal or elemental location for a tesnsor stress result.
# The default is the nodal location.

el_stress_result = solution.stress(location=post.locations.elemental)
nod_stress_result = solution.stress(location=post.locations.nodal)

###############################################################################
# Get information on tensor result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

el_stress = el_stress_result.tensor
nod_stress = nod_stress_result.tensor

el_field = el_stress[0]
el_field.location

nod_field = nod_stress[0]
nod_field.location

el_stress.get_max_data_at_field(0)
