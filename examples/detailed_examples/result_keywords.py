"""
.. _ref_result_keywords:

ANSYS DPF-POST Use result keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows how to get 
more precise result with keywords 
while calling a result from a 
solution object.
"""

###############################################################################
# Get started
from ansys.dpf import post


###############################################################################
# Get the solution object
import os
path = os.getcwd()
path += "/../../tests/testfiles/model_with_ns.rst"

solution = post.load_solution(path)

###############################################################################
# Get the keyword list
post.print_available_keywords()

###############################################################################
# Use a scoping on result
displacement_result = solution.displacement(location = post.locations.nodal, node_scoping = [1, 2, 3]) #default location is nodal
displacement = displacement_result.vector

displacement.get_data_at_field(0)

stress_with_elem_scop_result = solution.stress(location = post.locations.elemental_nodal, element_scoping = [1])
stress_with_elem_scop = stress_with_elem_scop_result.tensor

stress_with_elem_scop.get_data_at_field(0)

###############################################################################
# Use a named selection on result (MAPDL named_selection name must be written using UPPER CASE)
stress_on_ns_result = solution.stress(location = post.locations.elemental_nodal, named_selection="SELECTION") #MAPDL named_selection name must be written is upper case
stress_on_ns = stress_on_ns_result.tensor

stress_on_ns.num_fields

len(stress_on_ns[0])

###############################################################################
# Get a subresult
disp_x = displacement_result.x

stress_yz = stress_with_elem_scop_result.yz

stress_principal_1 = stress_on_ns_result.principal_3

###############################################################################
# Filter on a time / time_scoping / a set
print(solution.time_freq_support)

stress_on_time_1s_result = solution.stress(time = 1.)
stress_on_time_1s = stress_on_time_1s_result.tensor

displacement_on_set_1_result = solution.displacement(set = 1)
displacement_on_set_1 = displacement_on_set_1_result.vector

elastic_strain_with_time_scoping_result = solution.elastic_strain(time_scoping = [1, 3])
elastic_strain_with_time_scoping = elastic_strain_with_time_scoping_result.tensor

###############################################################################
# Make a grouping
displacement_result = solution.displacement(grouping = post.grouping.by_el_shape)
displacement_by_el_shape = displacement_result.vector

stress_result = solution.stress(mapdl_grouping = 186) #will filter only MAPDL elements of type solid 186
stress_on_solid_186 = stress_result.tensor

###############################################################################
# Manipulate the result object and change its definition
print(stress_on_ns_result)

stress_on_ns_result.definition.location

stress_on_ns_result.definition.location = post.locations.elemental
stress_on_ns_result.definition.time = 1.
stress_on_ns_elemental = stress_on_ns_result.tensor

print(stress_on_ns_result)

###############################################################################
# Use the misc. results
# Same keywords can be used here. For complex result, the keyword "phase" (float) can also be used.
stress_ratio = solution.misc.elemental_stress_ratio(node_scoping = [1, 32], time = 1.)
print(stress_ratio)

