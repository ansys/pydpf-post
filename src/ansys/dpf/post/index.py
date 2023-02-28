"""Module containing the ``Index`` class."""
from typing import List, Union

import ansys.dpf.core as dpf

location_to_label = {
    dpf.locations.nodal: "node",
    dpf.locations.elemental: "element",
    dpf.locations.elemental_nodal: "(element, node)",
    dpf.locations.overall: "overall",
    dpf.locations.time_freq_step: "step",
    dpf.locations.time_freq: "set",
}


class Index:
    """A Pandas style API to manipulate indexes."""

    def __init__(
        self,
        name: str,
        values: Union[List, None],
    ):
        """Creates an Index object to use in a DataFrame.

        Parameters
        ----------
        name:
            Name of the Index.
        values:
            Values taken by the Index.
        """
        self._name = name.replace(" ", "_")
        self._values = values
        self._dtype = None
        if values is not None:
            self._dtype = type(values[0])
        self._str = None

    def __repr__(self):
        """Representation of the Index."""
        return f'Index<"{self._name}", dtype={self._dtype}>'

    def __str__(self):
        """String representation of the Index."""
        return (
            f'Index "{self._name}" with '
            f"{len(self._values) if self._values is not None else 'uncounted'} "
            f"values of {self._dtype if self._dtype is not None else 'undetermined'} type."
        )

    @property
    def name(self):
        """Returns the name of the Index."""
        return self._name


class MultiIndex:
    """A Pandas style API to manipulate multi-indexes."""

    def __init__(
        self,
        indexes: List[Index],
    ):
        """Creates a MultiIndex from several Index objects.

        Parameters
        ----------
        indexes:
            List of class:`ansys.dpf.post.index.Index` objects.
        """
        self._indexes = indexes
        for i, index in enumerate(self._indexes):
            setattr(self, index.name, index)

    def __repr__(self):
        """Representation of the Index."""
        return "MultiIndex<".join([repr(index) + ", " for index in self._indexes]) + ">"

    def __str__(self):
        """String representation of the Index."""
        txt = f"MultiIndex with {len(self._indexes)} Index objects:\n"
        for index in self._indexes:
            txt += str(index) + "\n"
        return txt
