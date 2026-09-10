"""
.. _ref_introduction:

Introduction
============
This example shows how to load a solution and access some basic help.
"""
###############################################################################
# Perform required imports
# ------------------------
# Perform required imports.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Solution`` object
# -----------------------
# Get the ``Solution`` object. This example loads a simple file supplied with
# DPF-Post.

simple_bar = examples.simple_bar
solution = post.load_solution(simple_bar)

###############################################################################
# Use helpers
# -----------

###############################################################################
# Get solution information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Get solution information.

print(solution)

###############################################################################
# Get locations information
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Get location information. Locations are used to defined where
# the result must be set.

help(post.locations)

###############################################################################
# Get available keywords
# ~~~~~~~~~~~~~~~~~~~~~~
# Get the available keywords. Keywords can be used as attributes of a
# ``Result`` object.

post.print_available_keywords()
