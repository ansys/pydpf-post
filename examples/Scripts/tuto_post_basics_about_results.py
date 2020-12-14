from ansys.dpf import post

solution = post.load_solution("d:/rst/static.rst")

displacement_result = solution.displacement()
displacement = displacement_result.vector
displacement.num_fields
disp_data = displacement.get_data_at_field(0)
len(disp_data)
disp_data[1]
displacement.max_data
displacement.get_max_data_at_field(0)
displacement.min_data

el_stress_result = solution.stress(location = post.locations.elemental)
nod_stress_result = solution.stress(location = post.locations.nodal) #note: the default location is nodal
el_stress = el_stress_result.tensor
nod_stress = nod_stress_result.tensor
el_field = el_stress[0]
el_field.location
nod_field = nod_stress[0]
nod_field.location
el_stress.get_max_data_at_field(0)
