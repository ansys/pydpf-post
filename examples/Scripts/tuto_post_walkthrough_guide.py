from ansys.dpf import post

solution = post.load_solution("d:/rst/twobodies.rst")
mesh = solution.mesh
time_freq_support = solution.time_freq_support
post.print_available_keywords()

stress = solution.stress(location=post.locations.elemental, named_selection="SELECTION", time_scoping=[1]) 
print(stress)

stress.definition.location
stress.definition.named_selection
stress.definition.time_scoping

sx = stress.xx
sx.num_fields
sx_field = sx[0]
sx_data = sx.get_data_at_field(0)
len(sx_data)
sx.max
sx.max_data
sx.get_max_data_at_field(0)

s = stress.tensor
s_field = s[0]
s_data = sx.get_data_at_field(0)
len(s_data)
s.max
s.max_data
s.get_max_data_at_field(0)
s.min
s.min_data
s.get_min_data_at_field(0)
selection="SELECTION", time_scoping=[1]) 