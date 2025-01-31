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

"""Module containing helpers to build selections.

Selections
----------

"""
from ansys.dpf.post import selection


def select(
    time_freq_indexes=None,
    time_freq_sets=None,
    time_freq_values=None,
    named_selection_names=None,
    **kwargs,
):
    """Creates a ``Selection`` instance allowing to choose the domain on which to evaluate results.

    The results domain defines the time frequency and the spatial selection.

    Parameters
    ----------
    time_freq_indexes:
        Time/freq indexes to select.
    time_freq_sets:
        Time/freq sets to select.
    time_freq_values:
        Time/freq values to select.
    named_selection_names:
        Time/freq named selection to select.

    """
    current_selection = selection.Selection()
    if time_freq_indexes:
        current_selection.select_time_freq_indexes(time_freq_indexes)
    if time_freq_sets:
        current_selection.select_time_freq_sets(time_freq_sets)
    if time_freq_values:
        current_selection.select_time_freq_values(time_freq_values)
    if named_selection_names:
        current_selection.select_named_selection(named_selection_names)
    return current_selection
