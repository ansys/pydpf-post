from ansys.dpf import post

solution = post.load_solution("d:/rst/rth_complex.rst")
solution.has_complex_result()

disp_result = solution.displacement()
disp = disp_result.vector
disp.is_complex_result()
disp.num_fields

disp_at_phase = disp_result.vector_at_phase(39.)
disp_at_phase.max_data
disp_at_phase.num_fields
real_field = disp_result.vector_at_phase(0.)
img_field = disp_result.vector_at_phase(90.)

disp_ampl = disp_result.vector_amplitude
disp_ampl.max_data
disp_ampl.num_fields



