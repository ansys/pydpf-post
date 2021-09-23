"""
.. _ref_basics:

ANSYS DPF Basic Features of the POST API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows how to get a result object from ``ansys.dpf.post``
and its basic usage.
"""

###############################################################################
# **Get started**
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# **Get the solution object**: must be instantiated with the result filepath
example_path = examples.download_all_kinds_of_complexity()
solution = post.load_solution(example_path)

###############################################################################
# **Get a displacement result from the result object**
# First will be created a displacement result.

displacement_result = solution.displacement()
displacement = displacement_result.vector

###############################################################################
# **Use the displacement result**
displacement.num_fields
disp_data = displacement.get_data_at_field(0)
len(disp_data)

disp_data[1]

displacement.max_data
displacement.get_max_data_at_field(0)

displacement.min_data

###############################################################################
# **Get a stress result from the result object (nodal or elemental location)**
el_stress_result = solution.stress(location=post.locations.elemental)
nod_stress_result = solution.stress(
    location=post.locations.nodal
)  # note: the default location is nodal

el_stress = el_stress_result.tensor
nod_stress = nod_stress_result.tensor

el_field = el_stress[0]
el_field.location

nod_field = nod_stress[0]
nod_field.location

el_stress.get_max_data_at_field(0)
