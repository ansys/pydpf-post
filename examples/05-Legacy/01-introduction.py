# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
# PyDPF-Post.

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
