"""
.. _ref_introduction:

Load a solution
===============
This example shows how to use the legacy PyDPF-Post API to load a solution and
access some basic information.
"""
###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get by importing the DPF ``examples`` package.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Solution`` object
# -----------------------
# Get the ``Solution`` object. This example loads a simple file supplied with
# PyDPF-Post.

simple_bar = examples.simple_bar
solution = post.load_solution(simple_bar)

###############################################################################
# Use helpers
# -----------

###############################################################################
# Get solution information
# ~~~~~~~~~~~~~~~~~~~~~~~~

print(solution)

###############################################################################
# Get location information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Locations are used to defined where the result must be set.

help(post.locations)

###############################################################################
# Get available keywords
# ~~~~~~~~~~~~~~~~~~~~~~
# Keywords can be used as attributes of a ``Result`` object.

post.print_available_keywords()
