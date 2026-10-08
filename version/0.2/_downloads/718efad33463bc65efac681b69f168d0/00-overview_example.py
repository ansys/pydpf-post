"""
.. _ref_overview_example:

DPF-Post overview
=================
This example provides an overview of how you use DPF-Post.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. # This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Solution`` object
# -----------------------
# Get the ``Solution`` object that allows access to the result. The ``Solution``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

solution = post.load_solution(examples.multishells_rst)

###############################################################################
# Get mesh and time frequency support
# -----------------------------------
# Get the mesh and time frequency support. The mesh is the support of the model.
# The time frequency support is the time/frequency representation of the model.

mesh = solution.mesh
time_freq_support = solution.time_freq_support

###############################################################################
# Get ``Result`` object
# ---------------------
# Get a ``Result`` object from the ``Solution`` object. The ``Result`` object
# can be a stress, plastic strain, elastic strain, temperature, or displacement.

post.print_available_keywords()
stress = solution.stress(location=post.locations.nodal, time_scoping=[1])

# Both location and ``time_scoping`` are available in the definition of
# the stress result.

stress.definition.location
stress.definition.time_scoping

print(stress)

###############################################################################
# Compute data
# ------------
# Compute data.
#
# **SX subresult**
#
# This code gets the subresult ``SX``, which is the stress tensor in the XX direction.

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
#
# This code gets the minimum and maximum stresses at a field for all
# directions (``XX``, ``XY``, ``XZ``, ``XY``, ``YZ``, and ``XZ``.

s = stress.tensor
s_field = s[0]
s_data = sx.get_data_at_field(0)
print("Length of the data:", len(s_data))
print()
print("Maximum stress field:\n", s.max)
print()
print("Maximum data at stress field:", s.max_data)
print("Maximum stress tensors at field 0:\n", s.get_max_data_at_field(0))

print("Minimum stress field:\n", s.min)
print()
print("Minimum data at stress field:", s.min_data)
print("Minimum stress tensors at field 0:\n", s.get_min_data_at_field(0))


###############################################################################
# Plot result
# -----------
# Plot a result by using the ``plot_contour()`` method.

s.plot_contour()
