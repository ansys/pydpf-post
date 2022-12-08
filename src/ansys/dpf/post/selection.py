"""
Selection
---------

"""

from ansys.dpf.core import Scoping, time_freq_scoping_factory, Field, locations, operators, \
    Workflow, types
from ansys.dpf.core.field import _get_size_of_list
from typing import Union


class _WfNames:
    data_sources = "data_sources"
    scoping = "scoping"
    scoping_a = "scoping_a"
    scoping_b = "scoping_b"
    streams = "streams"


class Selection:
    """The ``Selection`` class helps defining the domain on which results are evaluated.
    The results domain defines the time frequency and the spatial selection.

    Parameters
    ----------


    """

    def __init__(self, server=None):
        self._server = server
        self._spatial_selection = SpatialSelection(server=server)
        self._time_freq_selection = TimeFreqSelection(server=server)
        self._qualifier_selection = None

    @property
    def time_freq_selection(self):
        """Returns the computed one based time scoping.

        Returns
        -------
        selection : TimeFreqSelection
        """
        return self._time_freq_selection

    @time_freq_selection.setter
    def time_freq_selection(self, value):
        self._time_freq_selection = value

    def select_time_freq_indices(self, time_freq_indices) -> None:
        """Select time frequency sets by their indices (0 based).

        Parameters
        ----------
        time_freq_indices: list[int]
        """
        self._time_freq_selection.select_time_freq_indices(time_freq_indices)

    def select_time_freq_sets(self, time_freq_sets) -> None:
        """Select time frequency sets by their cumulative sets (1 based).

        Parameters
        ----------
        time_freq_indices: list[int]
        """
        self._time_freq_selection.select_time_freq_sets(time_freq_sets)

    def select_time_freq_values(self, time_freq_values) -> None:
        """Select time frequency sets by their cumulative sets (1 based).

        Parameters
        ----------
        time_freq_indices: list[float], np.ndarray, ansys.dpf.core.Field
        """
        self._time_freq_selection.select_time_freq_values(time_freq_values)

    def select_named_selection(self, named_selection_name, location=None) -> None:
        """Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        named_selection_name : list[str]

        location: str, ansys.dpf.core.locations, optional
        """
        self._spatial_selection.select_named_selection(named_selection_name, location)

    def select_nodes(self, nodes: Union[list[int], Scoping]) -> None:
        """Select a mesh scoping with its node Ids.
        Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        nodes : list[int], Scoping
            node Ids.
        """
        self._spatial_selection.select_nodes(nodes)

    def select_elements(self, elements: Union[list[int], Scoping]) -> None:
        """Select a mesh scoping with its node Ids.
        Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        elements : list[int], Scoping
            element Ids.
        """
        self._spatial_selection.select_elements(elements)


class TimeFreqSelection:
    """Translate time frequency selection types (index, sets, time frequency values)
    into DPF known entities.
    """

    def __init__(self, server):
        self._server = server
        self._selection = Workflow(server=self._server)

    def select_time_freq_indices(self, time_freq_indices) -> None:
        """Select time frequency sets by their indices (0 based).

        Parameters
        ----------
        time_freq_indices: list[int]
        """
        time_freq_sets = [d + 1 for d in time_freq_indices]
        self.select_time_freq_sets(time_freq_sets)

    def select_time_freq_sets(self, time_freq_sets) -> None:
        """Select time frequency sets by their cumulative sets (1 based).

        Parameters
        ----------
        time_freq_indices: list[int]
        """
        sets = time_freq_scoping_factory.scoping_by_sets(time_freq_sets,
                                                         server=self._server)
        op = operators.utility.forward(
            sets,
            server=self._server
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_time_freq_values(self, time_freq_values) -> None:
        """Select time frequency sets by their cumulative sets (1 based).

        Parameters
        ----------
        time_freq_indices: list[float], np.ndarray, ansys.dpf.core.Field
        """
        if isinstance(time_freq_values, Field):
            values = time_freq_values
        else:
            time_freq_field = Field(location=locations.time_freq, server=self._server)
            time_freq_field.data = time_freq_values
            time_freq_field.scoping.ids = range(1, _get_size_of_list(time_freq_values))
            values = time_freq_field
        op = operators.utility.forward(
            values,
            server=self._server
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def evaluate_on(self, solution) -> Scoping:
        """Returns what is evaluated from the selections made on a given Solution.
        This scoping is internally used to evaluate result on the right spatial domain.

        Parameters
        ----------
        solution: Solution

        Returns
        -------
        Scoping
        """
        if self._selection is None:
            return None
        input_names = self._selection.input_names
        # TO DO: connect streams
        # if solution._model.metadata.streams_provider is not None and _WfNames.streams in input_names:
        #     self._selection.connect(_WfNames.streams,
        #                             solution._model.metadata.streams_provider.outputs.streams_container())
        # el
        if _WfNames.data_sources in input_names:
            self._selection.connect(_WfNames.data_sources, solution._model.metadata.data_sources)

        return self._selection.get_output(_WfNames.scoping, types.scoping)


class SpatialSelection:
    """Translate space selection types (node ids, elements ids, geometry selections, named selections...)
    into DPF known entities.
    """

    def __init__(self, selection=None, server=None):
        self._server = server
        self._selection = Workflow(server=self._server)

        if selection is not None:
            self.select_with_scoping(selection)

    def select_named_selection(self, named_selection_name, location=None) -> None:
        """Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        named_selection_names : list[str]

        location: str, ansys.dpf.core.locations, optional
        """
        op = operators.scoping.on_named_selection(
            requested_location=location,
            named_selection_name=named_selection_name,
            server=self._server
        )
        self._selection.add_operator(op)
        self._selection.set_input_name(_WfNames.data_sources, op.inputs.data_sources)
        self._selection.set_input_name(_WfNames.streams, op.inputs.streams_container)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.mesh_scoping)

    def select_with_scoping(self, scoping: Scoping):
        """Directly sets the scoping as the selection.

        Parameters
        ----------
        scoping: Scoping
        """
        if not isinstance(scoping, Scoping):
            raise TypeError(f"The input scoping is an instance of {str(type(scoping))} "
                            f"instead of an expected {str(Scoping)}.")

        op = operators.utility.forward(
            scoping,
            server=self._server
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_nodes(self, nodes: Union[list[int], Scoping]) -> None:
        """Select a mesh scoping with its node Ids.
        Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        nodes : list[int], Scoping
            node Ids.
        """
        if isinstance(nodes, Scoping):
            scoping = nodes
        else:
            scoping = Scoping(location=locations.nodal, ids=nodes, server=self._server)
        self.select_with_scoping(scoping)

    def select_elements(self, elements: Union[list[int], Scoping]) -> None:
        """Select a mesh scoping with its node Ids.
        Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        elements : list[int], Scoping
            element Ids.
        """
        if isinstance(elements, Scoping):
            scoping = elements
        else:
            scoping = Scoping(location=locations.elemental, ids=elements, server=self._server)
        self.select_with_scoping(scoping)

    def select_intersection_with(self, spatial_selection) -> None:
        """Keeps the Intersection between the current selection and the input one

        Parameters
        ----------
        spatial_selection: Selection, SpatialSelection, ansys.dpf.core.Scoping

        """
        if isinstance(spatial_selection, Selection):
            spatial_selection = spatial_selection._spatial_selection
        elif isinstance(spatial_selection, Scoping):
            spatial_selection = SpatialSelection(selection=spatial_selection, server=self._server)

        intersect_op = operators.scoping.intersect(server=self._server)

        new_wf = Workflow(self._server)
        new_wf.add_operator(intersect_op)
        new_wf.set_input_name(_WfNames.scoping_a, intersect_op.inputs.scopingA)
        new_wf.set_input_name(_WfNames.scoping_b, intersect_op.inputs.scopingB)
        new_wf.set_output_name(_WfNames.scoping, intersect_op.outputs.intersection)
        new_wf.connect_with(self._selection, {_WfNames.scoping: _WfNames.scoping_a})
        new_wf.connect_with(spatial_selection._selection, {_WfNames.scoping: _WfNames.scoping_b})
        self._selection = new_wf

    def evaluate_on(self, solution) -> Scoping:
        """Returns what is evaluated from the selections made on a given Solution.
        This scoping is internally used to evaluate result on the right spatial domain.

        Parameters
        ----------
        solution: Solution

        Returns
        -------
        Scoping
        """
        if self._selection is None:
            return None
        input_names = self._selection.input_names
        # TO DO: connect streams
        # if solution._model.metadata.streams_provider is not None and _WfNames.streams in input_names:
        #     self._selection.connect(_WfNames.streams,
        #                             solution._model.metadata.streams_provider.outputs.streams_container())
        # el
        if _WfNames.data_sources in input_names:
            self._selection.connect(_WfNames.data_sources, solution._model.metadata.data_sources)

        return self._selection.get_output(_WfNames.scoping, types.scoping)
