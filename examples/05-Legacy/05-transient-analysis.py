"""
.. _ref_trasient_analysis:

Transient analysis
==================
This example shows how you can post-process a result file for a transient analysis
using PyDPF-Post.
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
# Get the ``Solution`` object. This example loads a result file for a transient
# analysis computed in Ansys Mechanical.

solution = post.load_solution(examples.msup_transient)
print(solution)

###############################################################################
# Get ``Result`` objects
# ----------------------

###############################################################################
# Get displacement result
# ~~~~~~~~~~~~~~~~~~~~~~~
# Get the displacement ``Result`` object.

disp_result = solution.displacement()
disp = disp_result.vector

###############################################################################
# Check number of fields
# ~~~~~~~~~~~~~~~~~~~~~~
# Check the number of fields.

disp.num_fields

###############################################################################
# Get data from field
# ~~~~~~~~~~~~~~~~~~~
# Get data from a field.

disp.get_data_at_field(0)

###############################################################################
# Get maximum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the maximum data value over all fields.

disp.max_data

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the minimum data value over all fields.

disp.min_data

###############################################################################
# Get maximum data value over targeted field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the maximum data value over a targeted field.

disp.get_max_data_at_field(0)

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the minimum data value over all fields.

disp.get_min_data_at_field(0)

###############################################################################
# Get stress result
# -----------------
# Get the stress ``Result`` object for a tensor.

stress_result = solution.stress()
stress = stress_result.tensor

###############################################################################
# Check number of fields
# ~~~~~~~~~~~~~~~~~~~~~~
# Check the number of shell and solid elements in distinct fields.

stress.num_fields

###############################################################################
# Get shell field
# ~~~~~~~~~~~~~~~
# Get the shell field.

shell_field = stress[0]
shell_field.shell_layers

###############################################################################
# Get solid field
# ~~~~~~~~~~~~~~~
# Get the solid field.

solid_field = stress[0]

###############################################################################
# Plot contour
# ~~~~~~~~~~~~
# Plot the contour.

stress.plot_contour()

###############################################################################
# Get elastic strain result
# -------------------------
# Get an elastic strain result.

elastic_strain_result = solution.elastic_strain()
elastic_strain = elastic_strain_result.tensor
# shell and solid elements are in distinct fields.
elastic_strain.num_fields

###############################################################################
# If the result file contains results, you can use this method
# to get the elastic strain result.

print(solution.elastic_strain())
