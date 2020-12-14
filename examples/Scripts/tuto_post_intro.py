from ansys.dpf import post
solution = post.load_solution("d:/rst/static.rst")
print(solution)
print(solution.get_result_info())
help(post.locations)
post.print_available_keywords()