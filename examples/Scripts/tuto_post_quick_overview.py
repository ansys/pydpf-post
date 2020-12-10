from ansys.dpf import post

solution = post.load_solution("d:/rst/twobodies.rst")
mesh = solution.mesh
time_freq_support = solution.time_freq_support

stress = solution.stress(location=post.locations.elemental, time_scoping=[1]) 

sx = stress.xx
s2 = stress.principal_2
s = stress.tensor

sx.max_data
s2.num_fields
s.get_data_at_field(0) #equivalent to s[0].data

s.plot_contour("time", 1)

disp = solution.displacement()
d = disp.vector
d.plot_contour()

