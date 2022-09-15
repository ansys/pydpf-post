"""
.. _animations:

Animation of ResultData
=======================
This example shows how to create animations from a ResultData.
"""

from ansys.dpf import post
from ansys.dpf.post import examples

# Load the solution
solution = post.load_solution(examples.msup_transient)

# Get the ResultData for displacement
displacement_result = solution.displacement()
displacement = displacement_result.vector

# Create an animation using the animate method
displacement.animate()
