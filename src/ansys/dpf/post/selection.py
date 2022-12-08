"""
Selection
---------

"""
from __future__ import annotations
from ansys.dpf.post.solution import Solution
from ansys.dpf.core import Scoping, time_freq_scoping_factory, Field, locations, operators, \
    Workflow, types
from ansys.dpf.core.server_types import BaseServer
from ansys.dpf.core.server import get_or_create_server
from ansys.dpf.core.field import _get_size_of_list
from numpy import ndarray
from typing import Union, Optional


class _WfNames:
    data_sources = "data_sources"
    scoping = "scoping"
    scoping_a = "scoping_a"
    scoping_b = "scoping_b"
    streams = "streams"


class TimeFreqSelection:
    """Translate time/frequency selection types (index, sets, time/frequency values)
    into DPF known entities.
    """

    def __init__(self, server: Optional[BaseServer] = None):
        self._server = get_or_create_server(server)
        self._selection = Workflow(server=self._server)

    def select_time_freq_indices(self, time_freq_indices: list[int]) -> None:
        """Select time frequency sets by their indices (0 based).

        Parameters
        ----------
        time_freq_indices:
            Indices of the time-steps/frequencies to select.
        """
        time_freq_sets = [d + 1 for d in time_freq_indices]
        self.select_time_freq_sets(time_freq_sets)

    def select_time_freq_sets(self, time_freq_sets: list[int]) -> None:
        """Select time frequency sets by their cumulative sets (1 based).

        Parameters
        ----------
        time_freq_sets:
            Cumulative set indices of the time-steps/frequencies to select.
        """
        sets = time_freq_scoping_factory.scoping_by_sets(time_freq_sets,
                                                         server=self._server)
        op = operators.utility.forward(
            sets,
            server=self._server
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_time_freq_values(self,
                                time_freq_values: Union[list[float], ndarray, Field]
                                ) -> None:
        """Select time frequency sets by their values (1 based).

        Parameters
        ----------
        time_freq_values:
            Time/frequency values to select.
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

    def _evaluate_on(self, solution: Solution) -> Optional[Scoping]:
        """Returns what is evaluated from the selections made on a given Solution.
        This scoping is internally used to evaluate result on the right time/freq domain.

        Parameters
        ----------
        solution:
            DPF-Post Solution to evaluate the time/freq selection on.

        Returns
        -------
        Scoping:
            Resulting time/freq scoping. Returns ``None`` if no selection was made beforehand.
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

    def apply_to(self, solution: Solution) -> list[int]:
        """Performs the currently defined selection on the given Solution.

        Parameters
        ----------
        solution:
            PyDPF-Post Solution to apply the selection on.
        Returns
        -------
        IDs of the entities obtained after applying the selection.
        """
        scoping = self._evaluate_on(solution=solution)
        return scoping.ids


class SpatialSelection:
    """Translate space selection types (node ids, elements ids, geometry selections, named selections...)
    into DPF known entities.
    """

    def __init__(self, scoping: Optional[Scoping] = None,
                 server: Optional[BaseServer] = None):
        self._server = get_or_create_server(server)
        self._selection = Workflow(server=self._server)

        if scoping is not None:
            self.select_with_scoping(scoping)

    def select_named_selection(self,
                               named_selection: str,
                               location: Optional[Union[str, locations]] = None
                               ) -> None:
        """Select a mesh scoping corresponding to a named selection.

        Parameters
        ----------
        named_selection:
            Name of the named selection.

        location:
            Location to use (nodal, elemental...)
        """
        op = operators.scoping.on_named_selection(
            requested_location=location,
            named_selection_name=named_selection,
            server=self._server
        )
        self._selection.add_operator(op)
        self._selection.set_input_name(_WfNames.data_sources, op.inputs.data_sources)
        self._selection.set_input_name(_WfNames.streams, op.inputs.streams_container)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.mesh_scoping)

    def select_with_scoping(self, scoping: Scoping):
        """Directly sets the scoping as the spatial selection.

        Parameters
        ----------
        scoping:
            Scoping to use for spatial selection.
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
        """Select nodes using their IDs or a nodal mesh scoping.

        Parameters
        ----------
        nodes :
            node IDs or nodal mesh scoping.
        """
        if isinstance(nodes, Scoping):
            scoping = nodes
        else:
            scoping = Scoping(location=locations.nodal, ids=nodes, server=self._server)
        self.select_with_scoping(scoping)

    def select_elements(self, elements: Union[list[int], Scoping]) -> None:
        """Select elements using their IDs or an elemental mesh scoping.

        Parameters
        ----------
        elements :
            element IDs or elemental mesh scoping.
        """
        if isinstance(elements, Scoping):
            scoping = elements
        else:
            scoping = Scoping(location=locations.elemental, ids=elements, server=self._server)
        self.select_with_scoping(scoping)

    def intersect(self,
                  spatial_selection: Union[Selection, SpatialSelection, Scoping]
                  ) -> None:
        """Replaces the current selection by its intersection with the input.

        Parameters
        ----------
        spatial_selection:
            Spatial domain to intersect with.

        """
        if isinstance(spatial_selection, Selection):
            spatial_selection = spatial_selection.spatial_selection
        elif isinstance(spatial_selection, Scoping):
            spatial_selection = SpatialSelection(scoping=spatial_selection, server=self._server)

        intersect_op = operators.scoping.intersect(server=self._server)

        new_wf = Workflow(self._server)
        new_wf.add_operator(intersect_op)
        new_wf.set_input_name(_WfNames.scoping_a, intersect_op.inputs.scopingA)
        new_wf.set_input_name(_WfNames.scoping_b, intersect_op.inputs.scopingB)
        new_wf.set_output_name(_WfNames.scoping, intersect_op.outputs.intersection)
        new_wf.connect_with(self._selection, {_WfNames.scoping: _WfNames.scoping_a})
        new_wf.connect_with(spatial_selection._selection, {_WfNames.scoping: _WfNames.scoping_b})
        self._selection = new_wf

    def _evaluate_on(self, solution: Solution) -> Optional[Scoping]:
        """Performs the currently defined selection on the given Solution.
        This scoping is internally used to evaluate result on the right spatial domain.

        Parameters
        ----------
        solution:
            PyDPF-Post Solution to apply the selection on.
        Returns
        -------
        Scoping:
            Resulting time/freq scoping. Returns ``None`` if no selection was made beforehand.
        """
        if self._selection is None:
            return None
        input_names = self._selection.input_names
        if solution._model.metadata.streams_provider is not None and _WfNames.streams in input_names:
            self._selection.connect(_WfNames.streams,
                                    solution._model.metadata.streams_provider.outputs.streams_container())
        elif _WfNames.data_sources in input_names:
            self._selection.connect(_WfNames.data_sources, solution._model.metadata.data_sources)

        return self._selection.get_output(_WfNames.scoping, types.scoping)

    def apply_to(self, solution: Solution) -> list[int]:
        """Performs the currently defined selection on the given Solution.

        Parameters
        ----------
        solution:
            PyDPF-Post Solution to apply the selection on.
        Returns
        -------
        IDs of the entities obtained after applying the selection.
        """
        scoping = self._evaluate_on(solution=solution)
        return scoping.ids


class Selection:
    """The ``Selection`` class helps define the domain on which results are evaluated.
    The result domain defines the time/frequency and the spatial selection.

    Parameters
    ----------


    """

    def __init__(self, server: Optional[BaseServer] = None):
        self._server = get_or_create_server(server)
        self._spatial_selection = SpatialSelection(server=self._server)
        self._time_freq_selection = TimeFreqSelection(server=self._server)
        self._qualifier_selection = None

    @property
    def time_freq_selection(self) -> TimeFreqSelection:
        """Returns the computed one-based time scoping.

        Returns
        -------
        selection:
        """
        return self._time_freq_selection

    @time_freq_selection.setter
    def time_freq_selection(self, value: TimeFreqSelection):
        self._time_freq_selection = value

    @property
    def spatial_selection(self) -> SpatialSelection:
        """Returns the computed one-based mesh scoping.

        Returns
        -------
        selection:
        """
        return self._spatial_selection

    @spatial_selection.setter
    def spatial_selection(self, value: SpatialSelection):
        self._spatial_selection = value


    def select_time_freq_indices(self, time_freq_indices: list[int]) -> None:
        """Select time frequency sets by their indices (0 based).

        Parameters
        ----------
        time_freq_indices:
        """
        self._time_freq_selection.select_time_freq_indices(time_freq_indices)

    def select_time_freq_sets(self, time_freq_sets: list[int]) -> None:
        """Select time frequency sets by their cumulative sets (1 based).

        Parameters
        ----------
        time_freq_sets:
        """
        self._time_freq_selection.select_time_freq_sets(time_freq_sets)

    def select_time_freq_values(self,
                                time_freq_values: Union[]) -> None:
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
