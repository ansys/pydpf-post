"""
.. _ref_overview_example:

ANSYS DPF Overview of the POST Api
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows quickly how to use the POST Api.
It introduces the whole API. For more details, you 
can read the detailed tutorials.
"""

###############################################################################
# Start ansys.post module
# ~~~~~~~~~~~~~~~~~~~~~~~
# Here we load the post module.
from ansys.dpf import post

###############################################################################
# Prepare a solution object
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# We will first get the solution object that allows to access the result.
# It must be instantiated with the result filepath. 
import os
path = os.getcwd()
path += "/../../tests/testfiles/model_with_ns.rst"

solution = post.load_solution(path)

###############################################################################
# Get the mesh and the time_freq_support
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The mesh is the support of the model. 
# The time freq support is the time/frequence representation of the model.
mesh = solution.mesh
time_freq_support = solution.time_freq_support

###############################################################################
# Get a result object from the solution object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The result can be: stress, plastic_strain, elastic_strain, temperature, 
# displacement). 
post.print_available_keywords()
stress = solution.stress(location=post.locations.elemental, named_selection="SELECTION", time_scoping=[1]) 
print(stress)

stress.definition.location
stress.definition.named_selection
stress.definition.time_scoping

###############################################################################
# Compute data
# ~~~~~~~~~~~~
# Get a subresult (such as SX, or the stress tensor)

###############################################################################
# SX subresult
sx = stress.xx
sx.num_fields
sx_field = sx[0]
sx_data = sx.get_data_at_field(0)
len(sx_data)
sx.max
sx.max_data
sx.get_max_data_at_field(0)

###############################################################################
# Stress tensor result
s = stress.tensor
s_field = s[0]
s_data = sx.get_data_at_field(0)
len(s_data)
s.max
s.max_data
s.get_max_data_at_field(0)
s.min
s.min_data
s.get_min_data_at_field(0)


###############################################################################
# Display a result
# ~~~~~~~~~~~~~~~~
# Use <subresult>.plot_contour() to plot a result. 
s.plot_contour()




