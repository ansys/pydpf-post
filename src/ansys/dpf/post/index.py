"""Module containing the ``Index`` class."""
from typing import List, Union
import weakref

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
        values: Union[List, None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Creates an Index object to use in a DataFrame.

        Parameters
        ----------
        name:
            Name of the Index.
        values:
            Values taken by the Index.
        scoping:
            Scoping corresponding to this index to keep a weak reference.
        """
        self._name = name.replace(" ", "_")
        self._values = values
        self._dtype = None
        self._len = None
        self._scoping_ref = None
        # if scoping is None and values is None:
        #     raise ValueError("Arguments 'values' and 'scoping' cannot both be None.")
        if scoping is not None:
            self._scoping_ref = weakref.ref(scoping)
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
            f"values of {self._dtype if self._dtype is not None else 'undetermined'} type"
        )

    def __len__(self):
        """Returns the length of the index."""
        if self._len is not None:
            return self._len
        if self._scoping_ref is not None:
            return self._scoping_ref().size
        if self.values is not None:
            self._len = len(self.values)
            return self._len
        else:
            return None

    @property
    def name(self):
        """Returns the name of the Index."""
        return self._name

    @property
    def values(self):
        """Returns the values of the Index."""
        if self._values is None:
            self._evaluate_values()
        return self._values

    def _evaluate_values(self):
        """Evaluates the values of the Index."""
        if self._scoping_ref is not None:
            self._values = self._scoping_ref().ids


class ResultsIndex(Index):
    """Index class specific to results."""

    def __init__(
        self,
        values: List[str],
    ):
        """Initiate this class."""
        super().__init__(name="results", values=values, scoping=None)

    def __repr__(self):
        """Representation of the Index."""
        return f"ResultIndex<{self.values}>"


class MultiIndex:
    """A Pandas style API to manipulate multi-indexes."""

    def __init__(
        self,
        label_indexes: List[Index],
        results_index: Index,
    ):
        """Creates a MultiIndex from several Index objects.

        Parameters
        ----------
        label_indexes:
            List of class:`ansys.dpf.post.index.Index` objects relative to labels.
        results_index:
            Index relative to available results.
        """
        self._labels = label_indexes
        self._results = results_index
        for i, index in enumerate(self._labels):
            setattr(self, index.name, index)

    @property
    def labels(self):
        """Returns the list of Index objects for available labels."""
        return self._labels

    @property
    def results(self):
        """Returns the Index of available results."""
        return self._results

    def __repr__(self):
        """Representation of the Index."""
        return (
            "MultiIndex<".join([repr(index) + ", " for index in self.labels])
            + f"{repr(self.results)}>"
        )

    def __str__(self):
        """String representation of the Index."""
        txt = f"MultiIndex with {len(self.labels)} Label Index objects:\n"
        for index in self.labels:
            txt += str(index) + "\n"
        txt += f"and a ResultsIndex of size {len(self.results)}"
        return txt

    def __len__(self):
        """Returns the length of the MultiIndex (number of combinations)."""
        length = 1
        for index in self.labels:
            length = length * len(index)
        length = length * len(self.results)
        return length

    @property
    def label_names(self):
        """Returns a list with the name of each label Index."""
        return [index.name for index in self.labels]

    @property
    def result_names(self):
        """Returns a list with the available results."""
        return self.results.values
