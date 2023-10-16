"""Module containing the ``Index`` class and sub-classes.

Index
-----

"""
from abc import ABC
from typing import List, Union
import weakref

import ansys.dpf.core as dpf

from ansys.dpf import post


class ref_labels:
    """Reference naming for different common Indexes."""

    components = "components"
    results = "results"
    time = "time"
    modes = "modes"
    frequencies = "frequencies"
    set_ids = "set_ids"
    node_ids = "node_ids"
    face_ids = "face_ids"
    element_ids = "element_ids"
    cell_ids = "cell_ids"
    elemental_nodal = "element_ids"
    step = "step_ids"
    overall = "overall"
    element_node = "node"


location_to_label = {
    dpf.locations.nodal: ref_labels.node_ids,
    "cells": ref_labels.cell_ids,
    dpf.locations.elemental: ref_labels.element_ids,
    dpf.locations.elemental_nodal: ref_labels.elemental_nodal,
    dpf.locations.overall: ref_labels.overall,
    dpf.locations.time_freq_step: ref_labels.step,
    dpf.locations.time_freq: ref_labels.set_ids,
    None: "unknown",
}

# dpf.locations.faces is available only starting with ansys-dpg-gate 0.4.0
try:  # pragma: no cover
    location_to_label[dpf.locations.faces] = ref_labels.face_ids
except AttributeError as e:
    if "type object 'locations' has no attribute 'faces'" in str(e):
        pass
    else:
        raise e


class Index(ABC):
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
        if values is not None and len(values) > 0:
            self._dtype = type(values[0])
        self._str = None

    def __repr__(self):
        """Representation of the Index."""
        return f'{type(self).__name__}<name="{self._name}", dtype={self._dtype}>'

    def __str__(self):
        """String representation of the Index."""
        return (
            f'{type(self).__name__} "{self._name}" with '
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


class MeshIndex(Index):
    """Index class specific to mesh entities."""

    def __init__(
        self,
        location: Union[post.locations, str],
        values: Union[List[int], None] = None,
        scoping: Union[dpf.Scoping, None] = None,
        fc: Union[dpf.FieldsContainer, None] = None,
    ):
        """Initiate this class."""
        name = location_to_label[location]
        self.location = location
        if fc is None and values is None and scoping is None:
            raise ValueError(
                "Arguments 'values', 'scoping' and 'fc' cannot all be None."
            )
        if fc is not None:
            self._fc = weakref.ref(fc)
        super().__init__(name=name, values=values, scoping=scoping)
        self._dtype = int

    def _evaluate_values(self):
        """Evaluates the values of the MeshIndex."""
        if self._scoping_ref is not None:
            self._values = self._scoping_ref().ids
        else:
            # Merge the fields container scoping
            fc = self._fc()
            if fc is not None:
                merge_op = dpf.operators.utility.merge_scopings(server=fc._server)
                if float(fc._server.version) >= 5.0:
                    scopings = dpf.operators.utility.extract_scoping(
                        field_or_fields_container=fc,
                        server=fc._server,
                    ).outputs.mesh_scoping_as_scopings_container
                    merge_op.connect(0, scopings)
                else:
                    for i, f in enumerate(fc):
                        merge_op.connect(i, f.scoping)
                self._values = merge_op.eval().ids
            else:
                raise AttributeError(
                    "The FieldsContainer affiliated to the MeshIndex is no longer "
                    "available. Cannot evaluate Index.values."
                )


class ResultsIndex(Index):
    """Index class specific to results."""

    def __init__(
        self,
        values: List[str],
        units: Union[List[Union[str, None]], None] = None,
    ):
        """Initiate this class."""
        # Initialize results labels
        result_values = values
        if units is None:
            # If no units, initialize as list of None
            units = [None] * len(values)
        else:
            # Truncate units if longer than values
            if len(units) > len(values):
                units = units[: len(values)]
            # Update result labels with unit when necessary
            for i, unit in enumerate(units):
                result_values[i] += f" ({unit})" if unit is not None else ""
            # If units was shorter than values, extend with Nones
            if len(units) < len(values):
                units.extend([None] * (len(values) - len(units)))
        self._units = units
        super().__init__(name=ref_labels.results, values=result_values, scoping=None)

    @property
    def units(self) -> List[str]:
        """Return the list of units for this index."""
        return self._units

    def __repr__(self):
        """Representation of the Index."""
        return f"ResultIndex<{self.values}>"


class LabelIndex(Index):
    """Index class specific to labels."""

    def __init__(
        self,
        name: str,
        values: Union[List, None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Initiate this class."""
        super().__init__(name=name, values=values, scoping=scoping)

    def __repr__(self):
        """Representation of the Index."""
        return f"LabelIndex<name={self.name}, values={self.values}>"


class TimeIndex(Index):
    """Index class specific to time."""

    def __init__(
        self,
        values: Union[List, None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Initiate this class."""
        super().__init__(name=ref_labels.time, values=values, scoping=scoping)

    def __repr__(self):
        """Representation of the Index."""
        return f"TimeIndex<values={self.values}>"


class ModeIndex(Index):
    """Index class specific to modes."""

    def __init__(
        self,
        values: Union[List, None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Initiate this class."""
        super().__init__(name=ref_labels.modes, values=values, scoping=scoping)

    def __repr__(self):
        """Representation of the Index."""
        return f"ModeIndex<values={self.values}>"


class FrequencyIndex(Index):
    """Index class specific to frequency."""

    def __init__(
        self,
        values: Union[List, None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Initiate this class."""
        super().__init__(name=ref_labels.frequencies, values=values, scoping=scoping)

    def __repr__(self):
        """Representation of the Index."""
        return f"FrequencyIndex<values={self.values}>"


class SetIndex(LabelIndex):
    """Index class specific to set_ids."""

    def __init__(
        self,
        values: Union[List, None] = None,
        scoping: Union[dpf.Scoping, None] = None,
    ):
        """Initiate this class."""
        super().__init__(name=ref_labels.set_ids, values=values, scoping=scoping)

    def __repr__(self):
        """Representation of the Index."""
        return f"SetIndex<values={self.values}>"


class CompIndex(Index):
    """Index class specific to components."""

    def __init__(
        self,
        values: Union[List, None] = None,
    ):
        """Initiate this class."""
        super().__init__(name=ref_labels.components, values=values)


class ElementNodeIndex(Index):
    """Index class specific to elemental nodal results."""

    def __init__(
        self,
    ):
        """Initiate this class."""
        # We know there will be at least one node value per element.
        super().__init__(name=ref_labels.element_node, values=[0])


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
            Ordered list of class:`ansys.dpf.post.index.Index` objects.
        """
        self._indexes = indexes
        # self._labels = []
        # self._label_names = None
        # self._result_names = None
        for _, index in enumerate(self._indexes):
            setattr(self, index.name, index)

    @property
    def indexes(self) -> List[Index]:
        """Returns the list of Index in the MultiIndex."""
        return self._indexes

    # @property
    # def labels(self):
    #     """Returns the list of label Index objects."""
    #     return self._labels
    #
    # @property
    # def results(self):
    #     """Returns the Index of available results."""
    #     return self._results

    def __repr__(self):
        """Representation of the Index."""
        return f"MultiIndex<{self._indexes}>"

    # def __str__(self):
    #     """String representation of the Index."""
    #     txt = f"MultiIndex with {len(self)} Label Index objects:\n"
    #     for index in self._indexes:
    #         txt += str(index) + "\n"
    #     # txt += f"and a ResultsIndex of size {len(self.results)}"
    #     return txt

    def __len__(self):
        """Returns the number of Index objects in the MultiIndex."""
        return len(self._indexes)

    def __getitem__(self, item):
        """Get an Index in the MultiIndex."""
        return self._indexes[item]

    @property
    def names(self):
        """Returns a list with the name of each Index."""
        return [index.name for index in self._indexes]

    @property
    def results_index(self) -> Union[ResultsIndex, None]:
        """Returns the available ResultsIndex if present."""
        for index in self._indexes:
            if isinstance(index, ResultsIndex):
                return index
        return None

    @property
    def mesh_index(self) -> Union[MeshIndex, None]:
        """Returns the available ResultsIndex if present."""
        for index in self._indexes:
            if isinstance(index, MeshIndex):
                return index
        return None

    @property
    def set_index(self) -> Union[SetIndex, None]:
        """Returns the available SetIndex if present."""
        for index in self._indexes:
            if isinstance(index, SetIndex):
                return index
        return None

    # @property
    # def label_names(self):
    #     """Returns a list with the name of each label Index."""
    #     if self._label_names is None:
    #         self._label_names = [index.name for index in self._indexes]
    #     return
    #
    # @property
    # def result_names(self):
    #     """Returns a list with the available results."""
    #     return self.results.values
