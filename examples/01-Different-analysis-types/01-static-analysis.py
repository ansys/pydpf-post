"""
.. _ref_static_analysis:

Static Analysis
~~~~~~~~~~~~~~~
This tutorial shows how post-process a static analysis result file
using API of the POST module.
"""

###############################################################################
# Get started
# ~~~~~~~~~~~
from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get the solution object
# ~~~~~~~~~~~~~~~~~~~~~~~
# The following file is the result of a static analysis computed using
# Ansys Mechanical.
#
# Here we load the solution

example_path = examples.download_all_kinds_of_complexity()
solution = post.load_solution(example_path)
print(solution)

###############################################################################
# Get result objects
# ~~~~~~~~~~~~~~~~~~

###############################################################################
# Get a displacement result
# =========================

###############################################################################
# **Get the result**
disp_result = solution.displacement()
disp = disp_result.vector
print(disp)

###############################################################################
# **Check the number of fields**
disp.num_fields

###############################################################################
# **Get data from a field**
disp.get_data_at_field(0)

###############################################################################
# **Get maximum data value over all fields**
disp.max_data

###############################################################################
# **Get minimum data value over all fields**
disp.min_data

###############################################################################
# **Get maximum data value over a targeted field**
disp.get_max_data_at_field(0)

###############################################################################
# **Get minimum data value over all fields**
disp.get_min_data_at_field(0)

###############################################################################
# Get a stress result
# ===================

###############################################################################
# **Get the result**
stress_result = solution.stress()
stress = stress_result.tensor
# shell and solid elements are in distinct fields.
stress.num_fields

###############################################################################
# **Get the shell field**
shell_field = stress[0]
shell_field.shell_layers

###############################################################################
# **Get the solid field field**
solid_field = stress[1]

###############################################################################
# **Plot the contour**
stress.plot_contour()

###############################################################################
# Get an elastic_strain result
# ============================

###############################################################################
# **Get the result**
elastic_strain_result = solution.elastic_strain()
elastic_strain = elastic_strain_result.tensor
# shell and solid elements are in distinct fields.
elastic_strain.num_fields

###############################################################################
# **It is also possible to deal with plastic_strain and temperature this way.**
# The result file must contain those results.
#
# Here is the plastic strain result
print(solution.plastic_strain())

###############################################################################
# and the temperature result
print(solution.structural_temperature())
