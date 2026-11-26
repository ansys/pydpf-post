"""
.. _ref_result_on_path:

Request result on a specific path
=================================
This example shows how you can request a result on a
specific path of coordinates.
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
# Get the ``Solution`` object.

solution = post.load_solution(examples.static_rst)

###############################################################################
# Create coordinates array
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Create the coordinates array to request the result on.

coordinates = [[0.024, 0.03, 0.003]]
for i in range(1, 51):
    coord_copy = coordinates[0].copy()
    coord_copy[1] = coord_copy[0] + i * 0.001
    coordinates.append(coord_copy)

###############################################################################
# Create ``DpfPath`` object
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a ``DpfPath`` object.

path = post.create_path_on_coordinates(coordinates=coordinates)

###############################################################################
# Request result on path
# ~~~~~~~~~~~~~~~~~~~~~~
# Request the result on this path.

stress = solution.stress(path=path)

###############################################################################
# Plot result
# -----------
# Plot the result.

stress_eqv = stress.von_mises
stress_eqv.plot_contour()
