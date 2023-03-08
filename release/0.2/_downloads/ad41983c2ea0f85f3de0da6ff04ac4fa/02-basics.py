"""
.. _ref_basics:

Basic features
==============
This example shows you how you can get and use a ``Result`` object.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports.

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
# Get the displacement ``Result`` object.

displacement_result = solution.displacement()
displacement = displacement_result.vector

###############################################################################
# Get information on result
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Get information on the displacement result.

displacement.num_fields
disp_data = displacement.get_data_at_field(0)
len(disp_data)

disp_data[1]

displacement.max_data
displacement.get_max_data_at_field(0)

displacement.min_data

###############################################################################
# Get stress result
# ~~~~~~~~~~~~~~~~~
# Get the stress ``Result`` object for a tensor. You can get the nodal or
# elemental location. The default is the nodal location.

el_stress_result = solution.stress(location=post.locations.elemental)
nod_stress_result = solution.stress(location=post.locations.nodal)

###############################################################################
# Get information on result
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Get information on the stress result.

el_stress = el_stress_result.tensor
nod_stress = nod_stress_result.tensor

el_field = el_stress[0]
el_field.location

nod_field = nod_stress[0]
nod_field.location

el_stress.get_max_data_at_field(0)
