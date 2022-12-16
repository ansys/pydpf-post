"""Example of how to use the DPF-Post module.

Show how to load a solution and plot the displacement of nodes 4 and 6.
"""

from ansys.dpf.post import examples, load_solution

solution = load_solution(examples.static_rst)

disp = solution.displacement(nodes=[4, 6])

# plotting
disp.plot(cmap="viridis")

# numpy as array
disp.as_array()

# max
disp.max()

# As dataframe
disp.as_data_frame()

# Animates the displacement (this is static though)
disp.animate()
