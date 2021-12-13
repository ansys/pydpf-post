"""
.. _ref_result_on_path:

ANSYS DPF-Post request results on specific path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows how to request the result over a
specific path of coordinates.
"""

###############################################################################
# Get started

from ansys.dpf import post
from ansys.dpf.post import examples
solution = post.load_solution(examples.static_rst)

###############################################################################
# Create a coordinates array on which request the result
coordinates = [[0.024, 0.03, 0.003]]
for i in range(1, 51):
    coord_copy = coordinates[0].copy()
    coord_copy[1] = coord_copy[0] + i * 0.001
    coordinates.append(coord_copy)

###############################################################################
# Create a DpfPath object
path = post.create_path_on_coordinates(coordinates=coordinates)

###############################################################################
# Request the result on it
stress = solution.stress(path=path)

###############################################################################
# Plot the result
stress_eqv = stress.von_mises
stress_eqv.plot_contour()