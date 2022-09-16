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
print(solution)
print("==============================")
# Get the ResultData for displacement
print(post.print_available_keywords())
print("==============================")
displacement_result = solution.displacement()  # node_scoping=range(1, 394))
displacement = displacement_result.x

print(displacement)
fc = displacement.result_fields_container
print(fc)
field = fc[0]
print(field)
mesh = field.meshed_region
print(mesh)
print("====")
print(mesh.nodes)
print(mesh.nodes.scoping)
# print(mesh.nodes.scoping.ids)
# scoping = mesh.nodes.scoping
# print(scoping)
# print(scoping.ids)

# # Create an animation using the animate method
# displacement.animate()
#
# stress_result = solution.stress()
# stress = stress_result.von_mises
#
# stress.animate()
