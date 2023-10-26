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
