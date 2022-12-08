"""
tools
=====
This module holds factories to create geometry, selections...
These factories help creating Objects used to defined which results are evaluated.
"""

from ansys.dpf.post import selection


def select(time_freq_indexes=None, time_freq_sets=None, time_freq_values=None, named_selection_names=None, kwargs=None):
    """Creates a ``Selection`` instance allowing to choose the domain on which results are evaluated.
    The results domain defines the time frequency and the spatial selection.

    Parameters
    ----------

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
