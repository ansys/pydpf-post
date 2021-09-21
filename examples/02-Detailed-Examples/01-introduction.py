"""
.. _ref_introduction:

ANSYS DPF Introduction to POST API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial shows the first step to load the post module and to
access some basic help.
"""

###############################################################################
# **Get started**
from ansys.dpf.post import examples
from ansys.dpf import post

###############################################################################
# **Load a solution**: this must be instantiated with the result filepath.
# Load a simple  example file included with DPF-Post
simple_bar = examples.simple_bar
solution = post.load_solution(simple_bar)

###############################################################################
# Use the helpers
# ~~~~~~~~~~~~~~~

###############################################################################
# **Solution information**
print(solution)

###############################################################################
# **Locations information**: the locations are used to defined where
# the result must be set.
help(post.locations)

###############################################################################
# **Get the available keywords** that can be used as attribute of a
# **result object**.
post.print_available_keywords()
