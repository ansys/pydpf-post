"""
.. _ref_trasient_analysis:

Postprocess a result file for a transient analysis
==================================================
This example shows how to use the legacy PyDPF-Post API to postprocess a result
file for a transient analysis.
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

disp_result = solution.displacement()
disp = disp_result.vector

###############################################################################
# Get number of fields
# ~~~~~~~~~~~~~~~~~~~~~~

disp.num_fields

###############################################################################
# Get data from field
# ~~~~~~~~~~~~~~~~~~~

disp.get_data_at_field(0)

###############################################################################
# Get maximum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.max_data

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.min_data

###############################################################################
# Get maximum data value over targeted field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.get_max_data_at_field(0)

###############################################################################
# Get minimum data value over all fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp.get_min_data_at_field(0)

###############################################################################
# Get stress result for a tensor
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

shell_field = stress[0]
shell_field.shell_layers

###############################################################################
# Get solid field
# ~~~~~~~~~~~~~~~

solid_field = stress[0]

###############################################################################
# Plot contour
# ~~~~~~~~~~~~

stress.plot_contour()

###############################################################################
# Get elastic strain result
# -------------------------

elastic_strain_result = solution.elastic_strain()
elastic_strain = elastic_strain_result.tensor
# shell and solid elements are in distinct fields.
elastic_strain.num_fields

###############################################################################
# If the result file contains results, use this method
# to get the elastic strain result.

print(solution.elastic_strain())
