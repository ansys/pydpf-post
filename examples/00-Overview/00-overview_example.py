"""
.. _ref_overview_example:

ANSYS DPF Overview of the POST API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows quickly how to use the DPF-POST API.  For more
details, you can read the detailed tutorials in the examples gallery.
"""

###############################################################################
# Start ansys.post module
# =======================
# Here we load the post module.
from ansys.dpf import post

###############################################################################
# Prepare a solution object
# =========================
# We will first get the solution object that allows to access the
# result.  It must be instantiated with the result filepath like
# ``"C:/Users/user/my_result.rst"`` or on Linux ``"/home/user/my_result.rst"``
#
# In this example, we use a built-in example included with this
# module, which you can get by importing the examples module from
# ``ansys.dpf.post``
from ansys.dpf.post import examples

solution = post.load_solution(examples.multishells_rst)

###############################################################################
# Get the mesh and the time_freq_support
# ======================================
# The mesh is the support of the model.
# The time freq support is the time/frequence representation of the model.
mesh = solution.mesh
time_freq_support = solution.time_freq_support

###############################################################################
# Get a result object from the solution object
# ============================================
# The result can be: stress, plastic_strain, elastic_strain, temperature,
# displacement).
post.print_available_keywords()
stress = solution.stress(location=post.locations.nodal, time_scoping=[1])

# both the location and time_scoping are available in the definition of
# the stress result
stress.definition.location
stress.definition.time_scoping

print(stress)

###############################################################################
# Compute data
# ============
# Here is shown how to get a subresult such as SX, the stress tensor
# in the XX direction.
#
# **SX subresult**
sx = stress.xx
sx.num_fields
sx_field = sx[0]
sx_data = sx.get_data_at_field(0)
print("Length of the data:", len(sx_data))
print()
print("Maximum Stress Field:\n", sx.max)
print()
print("Maximum data at stress field:", sx.max_data)
print("Maximum SX at Field 0:", sx.get_max_data_at_field(0))

###############################################################################
# **Stress tensor result**
# Here we get the minimum and maximum stresses at a field for all
# directions, including 'XX', 'XY, 'XZ', 'XY', 'YZ', and 'XZ'.
s = stress.tensor
s_field = s[0]
s_data = sx.get_data_at_field(0)
print("Length of the data:", len(s_data))
print()
print("Maximum Stress Field:\n", s.max)
print()
print("Maximum data at stress field:", s.max_data)
print("Maximum Stress Tensors at Field 0:\n", s.get_max_data_at_field(0))

print("Minimum Stress Field:\n", s.min)
print()
print("Minimum data at stress field:", s.min_data)
print("Minimum Stress Tensors at Field 0:\n", s.get_min_data_at_field(0))


###############################################################################
# Display a result
# ================
# Use <subresult>.plot_contour() to plot a result.
s.plot_contour()
