"""Module containing selection classes.

Selection
---------

"""
from __future__ import annotations

from typing import TYPE_CHECKING, List

if TYPE_CHECKING:  # pragma: no cover
    from ansys.dpf.post.simulation import Simulation
    from ansys.dpf.post.mesh import Mesh

from typing import Union

from ansys.dpf.core import (
    Field,
    Scoping,
    Workflow,
    locations,
    operators,
    time_freq_scoping_factory,
    types,
)
from ansys.dpf.core.common import natures
from ansys.dpf.core.field import _get_size_of_list
from ansys.dpf.core.server import get_or_create_server
from ansys.dpf.core.server_types import BaseServer
from numpy import ndarray


class _WfNames:
    data_sources = "data_sources"
    scoping = "scoping"
    final_scoping = "final_scoping"
    scoping_a = "scoping_a"
    scoping_b = "scoping_b"
    streams = "streams"
    initial_mesh = "initial_mesh"
    mesh = "mesh"
    external_layer = "external_layer"
    skin = "skin"
    read_cyclic = "read_cyclic"
    cyclic_sectors_to_expand = "cyclic_sectors_to_expand"
    cyclic_phase = "cyclic_phase"
    result = "result"


def _is_model_cyclic(is_cyclic: str):
    return is_cyclic in ["single_stage", "multi_stage"]


class TimeFreqSelection:
    """Define a time/frequency selection.

    Translate time/frequency selection types (index, sets, time/frequency values) into DPF known
    entities.
    """

    def __init__(self, server: Union[BaseServer, None] = None):
        """Create a time/frequency selection.

        Parameters
        ----------
        server:
            DPF server on which to create the TimeFreqSelection.
        """
        self._server = get_or_create_server(server)
        self._selection = Workflow(server=self._server)

    def select_time_freq_indices(self, time_freq_indices: List[int]) -> None:
        """Select time frequency sets by their indices (zero-based indexing).

        Parameters
        ----------
        time_freq_indices:
            Indices of the time-steps/frequencies to select.
        """
        time_freq_sets = [d + 1 for d in time_freq_indices]
        self.select_time_freq_sets(time_freq_sets)

    def select_time_freq_sets(self, time_freq_sets: List[int]) -> None:
        """Select time frequency sets by their cumulative sets (one-based indexing).

        Parameters
        ----------
        time_freq_sets:
            Cumulative set indices of the time-steps/frequencies to select.
        """
        sets = time_freq_scoping_factory.scoping_by_sets(
            time_freq_sets, server=self._server
        )
        op = operators.utility.forward(sets, server=self._server)
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_load_steps(self, load_steps: List[int]) -> None:
        """Select a list of load steps (one-based indexing).

        Parameters
        ----------
        load_steps:
            IDs of the load steps to select.
        """
        sets = time_freq_scoping_factory.scoping_by_load_steps(load_steps=load_steps)
        op = operators.utility.forward(sets, server=self._server)
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_with_scoping(self, scoping: Scoping):
        """Directly sets the scoping as the time/freq selection.

        Parameters
        ----------
        scoping:
            Scoping to use for time/freq selection.
        """
        if not isinstance(scoping, Scoping):
            raise TypeError(
                f"The input scoping is an instance of {str(type(scoping))} "
                f"instead of an expected {str(Scoping)}."
            )

        op = operators.utility.forward(scoping, server=self._server)
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_time_freq_values(
        self, time_freq_values: Union[List[float], ndarray, Field]
    ) -> None:
        """Select time frequency sets by their values (one-based indexing).

        Parameters
        ----------
        time_freq_values:
            Time/frequency values to select.
        """
        if isinstance(time_freq_values, Field):
            values = time_freq_values
        else:
            time_freq_field = Field(
                location=locations.time_freq, server=self._server, nature=natures.scalar
            )
            time_freq_field.data = time_freq_values
            time_freq_field.scoping.ids = range(1, _get_size_of_list(time_freq_values))
            values = time_freq_field
        op = operators.utility.forward(values, server=self._server)
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def _evaluate_on(self, simulation: Simulation) -> Union[Scoping, None]:
        """Returns what is evaluated from the selections made on a given Simulation.

        This scoping is internally used to evaluate result on the right time/freq domain.

        Parameters
        ----------
        simulation:
            Simulation object to evaluate the time/freq selection on.

        Returns
        -------
        Scoping:
            Resulting time/freq scoping. Returns ``None`` if no selection was made beforehand.
        """
        if self._selection is None:
            return None
        input_names = self._selection.input_names
        # TO DO: connect streams
        # if solution._model.metadata.streams_provider is not None
        # and _WfNames.streams in input_names:
        #     self._selection.connect(_WfNames.streams,
        #                             solution._model.metadata.streams_provider.outputs.streams_container())
        # el
        if _WfNames.data_sources in input_names:
            self._selection.connect(
                _WfNames.data_sources, simulation._model.metadata.data_sources
            )

        return self._selection.get_output(_WfNames.scoping, types.scoping)

    def apply_to(self, simulation: Simulation) -> List[int]:
        """Performs the currently defined selection on the given Simulation.

        Parameters
        ----------
        simulation:
            PyDPF-Post Simulation to apply the selection on.

        Returns
        -------
        IDs of the entities obtained after applying the selection.
        """
        scoping = self._evaluate_on(simulation=simulation)
        return scoping.ids


class SpatialSelection:
    """Define a spatial selection.

    Translate space selection types (node ids, elements ids, geometry selections,
    named selections...) into DPF known entities.
    """

    def __init__(
        self,
        scoping: Union[Scoping, None] = None,
        server: Union[BaseServer, None] = None,
    ):
        """Create a SpatialSelection.

        Parameters
        ----------
        scoping:
            Scoping to use to define the SpatialSelection.
        server:
            DPF server on which to create the SpatialSelection.
        """
        self._server = get_or_create_server(server)
        self._selection = Workflow(server=self._server)

        if scoping is not None:
            self.select_with_scoping(scoping)

    def select_named_selection(
        self,
        named_selection: Union[str, List[str]],
        location: Union[str, locations, None] = None,
    ) -> None:
        """Select a mesh scoping corresponding to one or several named selections.

        Parameters
        ----------
        named_selection:
            Name of the named selection.

        location:
            Location of the mesh entities to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal` or
            `post.locations.elemental`.
        """
        if isinstance(named_selection, str):
            op = operators.scoping.on_named_selection(
                requested_location=location,
                named_selection_name=named_selection,
                server=self._server,
            )
            self._selection.add_operator(op)
            self._selection.set_input_name(
                _WfNames.data_sources, op.inputs.data_sources
            )
            # self._selection.set_input_name(
            #     _WfNames.streams, op.inputs.streams_container
            # )
            self._selection.set_output_name(_WfNames.scoping, op.outputs.mesh_scoping)
        else:
            op = operators.utility.merge_scopings(server=self._server)
            forward_ds = operators.utility.forward(any=None, server=self._server)
            forward_sc = operators.utility.forward(any=None, server=self._server)
            self._selection.add_operators([op, forward_ds, forward_sc])
            self._selection.set_input_name(_WfNames.data_sources, forward_ds.inputs.any)
            # self._selection.set_input_name(_WfNames.streams, forward_sc.inputs.any)
            for pin, ns in enumerate(named_selection):
                mesh_scoping_op = operators.scoping.on_named_selection(
                    requested_location=location,
                    named_selection_name=ns,
                    server=self._server,
                    # data_sources=forward_ds.outputs.any,
                    # streams_container=forward_sc.outputs.any,
                )
                self._selection.add_operator(mesh_scoping_op)
                mesh_scoping_op.connect(3, forward_sc.outputs.any)
                mesh_scoping_op.connect(4, forward_ds.outputs.any)
                op.connect(pin, mesh_scoping_op.outputs.mesh_scoping)
            self._selection.set_output_name(_WfNames.scoping, op.outputs.merged_scoping)

    def select_external_layer(
        self,
        location: Union[locations, str] = locations.elemental,
        result_native_location: Union[str, locations, None] = None,
        elements: Union[List[int], Scoping, None] = None,
        is_model_cyclic: str = "not_cyclic",
    ) -> None:
        """Select the external layer of the mesh.

        Select the external layer of the mesh (possibly on parts of the mesh scoped to the
        ``elements`` input). The mesh corresponding to this external layer is then used to
        display results and the nodes or elements of this submesh are used to scope result data
        extraction.

        Parameters
        ----------
        location:
            Location of the mesh entities to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal` or
            `post.locations.elemental`.
        result_native_location:
            Native (as found in the file) location of the output result. Used to pick
            the location of the scoping.
        elements:
            List of elements to use to compute the external layer,
            default is all the elements of the model.
        is_model_cyclic:
            Cyclic type: ``not_cyclic``, ``single_stage``, or ``multi_stage``. This value can
            be returned by the Operator ``operators.metadata.is_cyclic``. Used to get the skin
            on the expanded mesh.
        """
        if (
            result_native_location is not None
            and _is_model_cyclic(is_model_cyclic)
            and location != result_native_location
        ):
            location = result_native_location
        op = operators.mesh.external_layer(server=self._server)
        self._selection.add_operator(op)

        if _is_model_cyclic(is_model_cyclic):
            mesh_provider_cyc = operators.mesh.mesh_provider()
            self._selection.add_operator(mesh_provider_cyc)
            self._selection.set_input_name(_WfNames.read_cyclic, mesh_provider_cyc, 14)
            self._selection.set_input_name(
                _WfNames.cyclic_sectors_to_expand, mesh_provider_cyc, 18
            )
            self._selection.set_input_name(
                _WfNames.data_sources, mesh_provider_cyc.inputs.data_sources
            )
            self._selection.set_input_name(
                _WfNames.streams, mesh_provider_cyc.inputs.streams_container
            )
            if elements is not None:
                if not isinstance(elements, Scoping):
                    elements = Scoping(
                        server=self._server, ids=elements, location=locations.elemental
                    )
                mesh_by_scop_op = operators.mesh.from_scoping(
                    mesh=mesh_provider_cyc,
                    scoping=elements,
                    server=self._server,
                )
                self._selection.add_operator(mesh_by_scop_op)
                op.inputs.mesh.connect(mesh_by_scop_op)
            else:
                op.inputs.mesh.connect(mesh_provider_cyc)
            self._selection.set_input_name(
                _WfNames.initial_mesh, mesh_provider_cyc, 100
            )  # hack
        elif elements is not None:
            if not isinstance(elements, Scoping):
                elements = Scoping(
                    server=self._server, ids=elements, location=locations.elemental
                )
            mesh_by_scop_op = operators.mesh.from_scoping(
                scoping=elements, server=self._server
            )
            self._selection.add_operator(mesh_by_scop_op)
            self._selection.set_input_name(
                _WfNames.initial_mesh, mesh_by_scop_op.inputs.mesh
            )
            op.inputs.mesh.connect(mesh_by_scop_op)
        else:
            self._selection.set_input_name(_WfNames.initial_mesh, op.inputs.mesh)

        self._selection.set_output_name(_WfNames.mesh, op.outputs.mesh)
        self._selection.set_output_name(_WfNames.external_layer, op.outputs.mesh)
        if location == locations.nodal:
            self._selection.set_output_name(
                _WfNames.scoping, op.outputs.nodes_mesh_scoping
            )
        elif location == locations.elemental or location == locations.elemental_nodal:
            self._selection.set_output_name(
                _WfNames.scoping, op.outputs.elements_mesh_scoping
            )

    def select_skin(
        self,
        location: Union[locations, str] = locations.elemental,
        result_native_location: Union[str, locations, None] = None,
        elements: Union[List[int], Scoping, None] = None,
        is_model_cyclic: str = "not_cyclic",
    ) -> None:
        """Select the skin of the mesh.

        Select the skin of the mesh (possibly on parts of the mesh scoped to the
        ``elements`` input). The mesh corresponding to this skin is then used to
        display results and the nodes or elements of this submesh are used to scope result data
        extraction.

        Parameters
        ----------
        location:
            Location of the mesh entities to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal` or
            `post.locations.elemental`.
        result_native_location:
            Native (as found in the file) location of the output result. Used to pick
            the location of the scoping.
        elements:
            List of elements to use to compute the external layer,
            default is all the elements of the model. Getting the skin on a selection of
            elements for cyclic symmetry is not supported.
        is_model_cyclic:
            Cyclic type: ``not_cyclic``, ``single_stage``, or ``multi_stage``. This value can
            be returned by the Operator ``operators.metadata.is_cyclic``. Used to get the skin
            on the expanded mesh.
        """
        op = operators.mesh.skin(server=self._server)
        self._selection.add_operator(op)
        mesh_input = op.inputs.mesh

        if _is_model_cyclic(is_model_cyclic):
            mesh_provider_cyc = operators.mesh.mesh_provider()
            self._selection.add_operator(mesh_provider_cyc)
            self._selection.set_input_name(_WfNames.read_cyclic, mesh_provider_cyc, 14)
            self._selection.set_input_name(
                _WfNames.cyclic_sectors_to_expand, mesh_provider_cyc, 18
            )
            self._selection.set_input_name(
                _WfNames.data_sources, mesh_provider_cyc.inputs.data_sources
            )
            self._selection.set_input_name(
                _WfNames.streams, mesh_provider_cyc.inputs.streams_container
            )
            if elements is not None:
                raise ValueError(
                    "Getting the skin an a selection of elements is not "
                    "supported for cyclic symmetry"
                )
                if not isinstance(elements, Scoping):
                    elements = Scoping(
                        server=self._server, ids=elements, location=locations.elemental
                    )
                mesh_by_scop_op = operators.mesh.from_scoping(
                    mesh=mesh_provider_cyc,
                    scoping=elements,
                    server=self._server,
                )
                self._selection.add_operator(mesh_by_scop_op)
                op.inputs.mesh.connect(mesh_by_scop_op)
            else:
                op.inputs.mesh.connect(mesh_provider_cyc)
            self._selection.set_input_name(
                _WfNames.initial_mesh, mesh_provider_cyc, 100
            )  # hack
            mesh_input = None

        elif elements is not None:
            if not isinstance(elements, Scoping):
                elements = Scoping(
                    server=self._server, ids=elements, location=locations.elemental
                )
            mesh_by_scop_op = operators.mesh.from_scoping(
                scoping=elements, server=self._server
            )
            self._selection.add_operator(mesh_by_scop_op)
            mesh_input = mesh_by_scop_op.inputs.mesh
            op.inputs.mesh.connect(mesh_by_scop_op)

        if mesh_input is not None:
            self._selection.set_input_name(_WfNames.initial_mesh, mesh_input)
            if location == result_native_location:
                self._selection.set_output_name(_WfNames.mesh, op.outputs.mesh)
        self._selection.set_output_name(_WfNames.skin, op.outputs.mesh)
        if location == locations.nodal and result_native_location == locations.nodal:
            self._selection.set_output_name(
                _WfNames.scoping, op.outputs.nodes_mesh_scoping
            )

        elif not _is_model_cyclic(is_model_cyclic) and (
            result_native_location == locations.elemental
            or result_native_location == locations.elemental_nodal
        ):
            transpose_op = operators.scoping.transpose(
                mesh_scoping=op.outputs.nodes_mesh_scoping, server=self._server
            )
            self._selection.add_operator(transpose_op)
            self._selection.set_input_name(
                _WfNames.initial_mesh, transpose_op.inputs.meshed_region
            )
            self._selection.set_output_name(
                _WfNames.scoping, transpose_op.outputs.mesh_scoping_as_scoping
            )

    def select_with_scoping(self, scoping: Scoping):
        """Directly sets the scoping as the spatial selection.

        Parameters
        ----------
        scoping:
            Scoping to use for spatial selection.
        """
        if not isinstance(scoping, Scoping):
            raise TypeError(
                f"The input scoping is an instance of {str(type(scoping))} "
                f"instead of an expected {str(Scoping)}."
            )

        op = operators.utility.forward(scoping, server=self._server)
        self._selection.add_operator(op)
        self._selection.set_output_name(_WfNames.scoping, op.outputs.any)

    def select_nodes(self, nodes: Union[List[int], Scoping]) -> None:
        """Select nodes using their IDs or a nodal mesh scoping.

        Parameters
        ----------
        nodes :
            node IDs or nodal mesh scoping.
        """
        if isinstance(nodes, Scoping):
            scoping = nodes
        else:
            from ansys.dpf.core.dpf_array import DPFArray

            if isinstance(nodes, DPFArray):
                nodes = nodes.tolist()
            scoping = Scoping(location=locations.nodal, ids=nodes, server=self._server)
        self.select_with_scoping(scoping)

    def select_nodes_of_elements(
        self,
        elements: Union[List[int], Scoping],
        mesh: Mesh,
    ) -> None:
        """Select all nodes of elements using the elements' IDs or an elemental mesh scoping.

        Parameters
        ----------
        elements:
            element IDs or elemental mesh scoping.
        mesh:
            Mesh containing the necessary connectivity.
        """
        if isinstance(elements, Scoping):
            scoping = elements
        else:
            scoping = Scoping(
                location=locations.elemental, ids=elements, server=self._server
            )

        op = operators.scoping.transpose(
            mesh_scoping=scoping,
            meshed_region=mesh._meshed_region,
            inclusive=0,
            requested_location=locations.nodal,
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(
            _WfNames.scoping, op.outputs.mesh_scoping_as_scoping
        )

    def select_nodes_of_faces(
        self,
        faces: Union[List[int], Scoping],
        mesh: Mesh,
    ) -> None:
        """Select all nodes of faces using the faces' IDs or a faces mesh scoping.

        Parameters
        ----------
        faces:
            face IDs or faces mesh scoping.
        mesh:
            Mesh containing the necessary connectivity.
        """
        if isinstance(faces, Scoping):
            scoping = faces
        else:
            scoping = Scoping(location=locations.faces, ids=faces, server=self._server)

        op = operators.scoping.transpose(
            mesh_scoping=scoping,
            meshed_region=mesh._meshed_region,
            inclusive=0,
            requested_location=locations.nodal,
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(
            _WfNames.scoping, op.outputs.mesh_scoping_as_scoping
        )

    def select_faces_of_elements(
        self,
        elements: Union[List[int], Scoping],
        mesh: Mesh,
    ) -> None:
        """Select all faces of elements using the elements' IDs or an elemental mesh scoping.

        Parameters
        ----------
        elements:
            element IDs or elemental mesh scoping.
        mesh:
            Mesh containing the necessary connectivity.
        """
        if isinstance(elements, Scoping):
            scoping = elements
        else:
            scoping = Scoping(
                location=locations.elemental, ids=elements, server=self._server
            )

        op = operators.scoping.transpose(
            mesh_scoping=scoping,
            meshed_region=mesh._meshed_region,
            inclusive=0,
            requested_location=locations.faces,
        )
        self._selection.add_operator(op)
        self._selection.set_output_name(
            _WfNames.scoping, op.outputs.mesh_scoping_as_scoping
        )

    def select_faces(self, faces: Union[List[int], Scoping]) -> None:
        """Select faces using their IDs or a faces mesh scoping.

        Parameters
        ----------
        faces :
            face IDs or faces mesh scoping.
        """
        if isinstance(faces, Scoping):
            scoping = faces
        else:
            scoping = Scoping(location=locations.faces, ids=faces, server=self._server)
        self.select_with_scoping(scoping)

    def select_elements(self, elements: Union[List[int], Scoping]) -> None:
        """Select elements using their IDs or an elemental mesh scoping.

        Parameters
        ----------
        elements :
            element IDs or elemental mesh scoping.
        """
        if isinstance(elements, Scoping):
            scoping = elements
        else:
            scoping = Scoping(
                location=locations.elemental, ids=elements, server=self._server
            )
        self.select_with_scoping(scoping)

    def intersect(
        self, spatial_selection: Union[Selection, SpatialSelection, Scoping]
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
            spatial_selection = SpatialSelection(
                scoping=spatial_selection, server=self._server
            )

        intersect_op = operators.scoping.intersect(server=self._server)

        new_wf = Workflow(self._server)
        new_wf.add_operator(intersect_op)
        new_wf.set_input_name(_WfNames.scoping_a, intersect_op.inputs.scopingA)
        new_wf.set_input_name(_WfNames.scoping_b, intersect_op.inputs.scopingB)
        new_wf.set_output_name(_WfNames.scoping, intersect_op.outputs.intersection)
        new_wf.connect_with(self._selection, {_WfNames.scoping: _WfNames.scoping_a})
        new_wf.connect_with(
            spatial_selection._selection, {_WfNames.scoping: _WfNames.scoping_b}
        )
        self._selection = new_wf

    def _evaluate_on(self, simulation: Simulation) -> Union[Scoping, None]:
        """Performs the currently defined selection on the given Simulation.

        This scoping is internally used to evaluate result on the right spatial domain.

        Parameters
        ----------
        simulation:
            PyDPF-Post Simulation to apply the selection on.

        Returns
        -------
        Scoping:
            Resulting time/freq scoping. Returns ``None`` if no selection was made beforehand.

        """
        if self._selection is None:
            return None
        input_names = self._selection.input_names
        if (
            simulation._model.metadata.streams_provider is not None
            and _WfNames.streams in input_names
        ):
            self._selection.connect(
                _WfNames.streams,
                simulation._model.metadata.streams_provider.outputs.streams_container(),
            )
        elif _WfNames.data_sources in input_names:
            self._selection.connect(
                _WfNames.data_sources, simulation._model.metadata.data_sources
            )

        return self._selection.get_output(_WfNames.scoping, types.scoping)

    def apply_to(self, simulation: Simulation) -> List[int]:
        """Performs the currently defined selection on the given Simulation.

        Parameters
        ----------
        simulation:
            PyDPF-Post Simulation to apply the selection on.

        Returns
        -------
        IDs of the entities obtained after applying the selection.

        """
        scoping = self._evaluate_on(simulation=simulation)
        return scoping.ids

    @property
    def requires_mesh(self) -> bool:
        """Whether the selection workflow requires a ``mesh`` as an input or not."""
        return _WfNames.initial_mesh in self._selection.input_names

    @property
    def outputs_mesh(self) -> bool:
        """Whether the selection workflow as an output named ``mesh``."""
        return _WfNames.mesh in self._selection.output_names

    def requires_manual_averaging(
        self,
        location: Union[str, locations],
        result_native_location: Union[str, locations],
        is_model_cyclic: str = "not_cyclic",
    ) -> bool:
        """Whether the selection workflow requires to manually build the averaging workflow."""
        output_names = self._selection.output_names
        is_model_cyclic = is_model_cyclic in ["single_stage", "multi_stage"]
        if (
            _WfNames.external_layer in output_names
            and is_model_cyclic
            and location != result_native_location
        ):
            return True
        elif _WfNames.skin in output_names and (
            result_native_location == locations.elemental
            or result_native_location == locations.elemental_nodal
        ):
            return True
        return False


class Selection:
    """The ``Selection`` class helps define the domain on which results are evaluated.

    The result domain defines the time/frequency and the spatial selection.
    """

    def __init__(self, server: Union[BaseServer, None] = None):
        """Instantiate a Selection.

        Parameters
        ----------
        server:
            DPF server on which to instantiate.
        """
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
            Time/freq selection of the selection.
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
            Spatial selection of the selection.
        """
        return self._spatial_selection

    @spatial_selection.setter
    def spatial_selection(self, value: SpatialSelection):
        self._spatial_selection = value

    def select_time_freq_indices(self, time_freq_indices: List[int]) -> None:
        """Select time frequency sets by their indices (zero-based indexing).

        Parameters
        ----------
        time_freq_indices:
            Time/freq indexes to select.
        """
        self._time_freq_selection.select_time_freq_indices(time_freq_indices)

    def select_time_freq_sets(self, time_freq_sets: Union[List[int], int]) -> None:
        """Select time frequency sets by their cumulative sets (one-based indexing).

        Parameters
        ----------
        time_freq_sets:
            Time/freq sets to select.
        """
        if isinstance(time_freq_sets, int):
            time_freq_sets = [time_freq_sets]
        self._time_freq_selection.select_time_freq_sets(time_freq_sets)

    def select_time_freq_values(
        self, time_freq_values: Union[List[float], ndarray, Field]
    ) -> None:
        """Select time frequency sets by their values.

        Parameters
        ----------
        time_freq_values:
            Time/freq values to select.
        """
        self._time_freq_selection.select_time_freq_values(time_freq_values)

    def select_named_selection(
        self,
        named_selection: Union[str, List[str]],
        location: Union[str, locations, None] = None,
    ) -> None:
        """Select a mesh scoping corresponding to one or several named selections.

        Parameters
        ----------
        named_selection:
            Named selection to select.
        location:
            Location of the mesh entities to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal` or
            `post.locations.elemental`.
        """
        self._spatial_selection.select_named_selection(named_selection, location)

    def select_nodes(self, nodes: Union[List[int], Scoping]) -> None:
        """Select a mesh scoping with its node IDs.

        Select a mesh scoping corresponding to a list of node IDs.

        Parameters
        ----------
        nodes:
            node IDs.
        """
        self._spatial_selection.select_nodes(nodes)

    def select_faces(self, faces: Union[List[int], Scoping]) -> None:
        """Select a mesh scoping with its face IDs.

        Select a mesh scoping corresponding to a list of face IDs.

        Parameters
        ----------
        faces:
            face IDs.
        """
        self._spatial_selection.select_faces(faces)

    def select_elements(self, elements: Union[List[int], Scoping]) -> None:
        """Select a mesh scoping with its element IDs.

        Select a mesh scoping corresponding to a list of element IDs.

        Parameters
        ----------
        elements:
            element IDs.
        """
        self._spatial_selection.select_elements(elements)

    def select_nodes_of_elements(
        self, elements: Union[List[int], Scoping], mesh: Mesh
    ) -> None:
        """Select nodes belonging to elements defined by their IDs.

        Select a nodal mesh scoping corresponding to elements.

        Parameters
        ----------
        elements:
            element IDs.
        mesh:
            Mesh containing the connectivity.
        """
        self._spatial_selection.select_nodes_of_elements(elements, mesh)

    def select_nodes_of_faces(
        self, faces: Union[List[int], Scoping], mesh: Mesh
    ) -> None:
        """Select nodes belonging to faces defined by their IDs.

        Select a nodal mesh scoping corresponding to faces.

        Parameters
        ----------
        faces:
            face IDs.
        mesh:
            Mesh containing the connectivity.
        """
        self._spatial_selection.select_nodes_of_faces(faces, mesh)

    def select_faces_of_elements(
        self, elements: Union[List[int], Scoping], mesh: Mesh
    ) -> None:
        """Select faces belonging to elements defined by their IDs.

        Select a faces mesh scoping corresponding to elements.

        Parameters
        ----------
        elements:
            element IDs.
        mesh:
            Mesh containing the connectivity.
        """
        self._spatial_selection.select_faces_of_elements(elements, mesh)

    def select_external_layer(
        self,
        location: Union[locations, str] = locations.elemental,
        result_native_location: Union[str, locations, None] = None,
        elements: Union[List[int], Scoping, None] = None,
        is_model_cyclic: str = "not_cyclic",
    ) -> None:
        """Select the external layer of the mesh.

        Select the external layer of the mesh (possibly on parts of the mesh scoped to the
        ``elements`` input). The mesh corresponding to this external layer are then used to
        display results and the nodes or elements of this submesh are used to scope result data
        extraction.

        Parameters
        ----------
        location:
            Location of the mesh entities to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal` or
            `post.locations.elemental`.
        result_native_location:
            Native (as found in the file) location of the output result. Used to pick
            the location of the scoping.
        elements:
            List of elements to use to compute the external layer,
            default is all the elements of the model.
        is_model_cyclic:
            Cyclic type: ``not_cyclic``, ``single_stage``, or ``multi_stage``. This value can
            be returned by the Operator ``operators.metadata.is_cyclic``. Used to get the skin
            on the expanded mesh.
        """
        self._spatial_selection.select_external_layer(
            elements=elements,
            location=location,
            result_native_location=result_native_location,
            is_model_cyclic=is_model_cyclic,
        )

    def select_skin(
        self,
        location: Union[locations, str] = locations.elemental,
        result_native_location: Union[str, locations, None] = None,
        elements: Union[List[int], Scoping, None] = None,
        is_model_cyclic: str = "not_cyclic",
    ) -> None:
        """Select the skin of the mesh.

        Select the skin of the mesh (possibly on parts of the mesh scoped to the
        ``elements`` input). The mesh corresponding to this skin are then used to
        display results and the nodes or elements of this submesh are used to scope result data
        extraction.

        Parameters
        ----------
        location:
            Location of the mesh entities to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal` or
            `post.locations.elemental`.
        result_native_location:
            Native (as found in the file) location of the output result. Used to pick
            the location of the scoping.
        elements:
            List of elements to use to compute the external layer,
            default is all the elements of the model. Getting the skin on a selection of
            elements for cyclic symmetry is not supported.
        is_model_cyclic:
            Cyclic type: ``not_cyclic``, ``single_stage``, or ``multi_stage``. This value can
            be returned by the Operator ``operators.metadata.is_cyclic``. Used to get the skin
            on the expanded mesh.
        """
        self._spatial_selection.select_skin(
            elements=elements,
            result_native_location=result_native_location,
            location=location,
            is_model_cyclic=is_model_cyclic,
        )

    @property
    def requires_mesh(self) -> bool:
        """Whether the selection workflow requires a ``initial_mesh`` as an input or not."""
        return self._spatial_selection.requires_mesh

    @property
    def outputs_mesh(self) -> bool:
        """Whether the selection workflow as an output named ``mesh``."""
        return self._spatial_selection.outputs_mesh

    def requires_manual_averaging(
        self,
        location: Union[str, locations],
        result_native_location: Union[str, locations],
        is_model_cyclic: str = "not_cyclic",
    ) -> bool:
        """Whether the selection workflow requires to manually build the averaging workflow."""
        return self._spatial_selection.requires_manual_averaging(
            location=location,
            result_native_location=result_native_location,
            is_model_cyclic=is_model_cyclic,
        )
