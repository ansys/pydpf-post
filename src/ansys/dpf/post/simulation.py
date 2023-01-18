"""Module containing the ``Simulation`` class."""
from abc import ABC
import re
from typing import List, Union

from ansys.dpf.core import DataSources, Model
from ansys.dpf.core.plotter import DpfPlotter

from ansys.dpf import core
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.selection import Selection


class Simulation(ABC):
    """Base class of all PyDPF-Post simulation types."""

    _component_id_to_str = {
        "1": "X",
        "2": "Y",
        "3": "Z",
        "4": "XY",
        "5": "YZ",
        "6": "XZ",
        "X": "X",
        "Y": "Y",
        "Z": "Z",
        "XY": "XY",
        "YZ": "YZ",
        "XZ": "XZ",
    }

    def __init__(self, data_sources: DataSources, model: Model):
        """Initialize the simulation using a ``dpf.core.Model`` object."""
        self._model = model
        self._data_sources = data_sources
        self._geometries = []
        self._boundary_conditions = []
        self._loads = []
        self._active_selection = None
        self._named_selections = None
        self._mesh = None

    @property
    def results(self) -> List[str]:
        r"""Available results.

        Returns a list of available results as strings.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.results) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['displacement\nOperator name: "U"\n...Units: degc\n']
        """
        return [
            str(result) for result in self._model.metadata.result_info.available_results
        ]

    @property
    def geometries(self):
        """List of constructed geometries in the simulation.

        Returns a list of geometry objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.geometries) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._geometries

    @property
    def boundary_conditions(self):
        """List of boundary conditions in the simulation.

        Returns a list of boundary_condition objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.boundary_conditions) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._boundary_conditions

    @property
    def loads(self):
        """List of loads in the simulation.

        Returns a list of load objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.loads) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._loads

    @property
    def mesh(self) -> Mesh:
        """Mesh representation of the model.

        Returns a :class:`ansys.dpf.post.mesh.Mesh` object.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.mesh.Mesh object at ...>
        """
        if self._mesh is None:
            self._mesh = Mesh(self._model.metadata.meshed_region)
        return self._mesh

    @property
    def named_selections(self) -> List[str]:
        """List of named selections in the simulation.

        Returns a list of named selections names.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.named_selections) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['_FIXEDSU']
        """
        if self._named_selections is None:
            self._named_selections = self._model.metadata.available_named_selections
        return self._named_selections

    def plot(
        self,
        mesh: bool = True,
        geometry: bool = True,
        loads: bool = True,
        boundary_conditions: bool = True,
    ):
        """General plot of the simulation object.

        Plots by default the complete mesh contained in the simulation,
        as well as a representation of the constructed geometry,
        the loads, and the boundary conditions currently defined.
        Each representation can be deactivated with its respective boolean argument.

        Args:
            mesh:
                Whether to plot the mesh representation.
            geometry:
                Whether to plot the geometries.
            loads:
                Whether to plot the loads.
            boundary_conditions:
                Whether to plot the boundary conditions.

        Returns
        -------
            Returns a plotter instance of the active visualization backend.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> simulation.plot() # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        """
        plt = DpfPlotter()
        if mesh:
            plt.add_mesh(self.mesh._meshed_region)
        if geometry:
            for geom in self.geometries:
                getattr(plt, "add_" + str(type(geom).__name__).lower())(geom)
        if loads:
            pass
        if boundary_conditions:
            pass
        plt.show_figure()

    @property
    def active_selection(self) -> Selection:
        """Active selection used by default for result queries.

        Returns a :object:`ansys.dpf.post.selection.Selection` object.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.activate_selection(selection=selection)
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        """
        return self._active_selection

    def activate_selection(self, selection: Selection):
        """Sets a selection as active on the simulation.

        Activating a given selection on a simulation means it is used
        as a default selection/filter in further result queries.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.activate_selection(selection=selection)
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        """
        self._active_selection = selection

    def deactivate_selection(self):
        """Deactivate the currently active selection.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.activate_selection(selection=selection)
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        >>> simulation.deactivate_selection()
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        None
        """
        self._active_selection = None

    @property
    def _time_frequencies(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._model.metadata.time_freq_support

    @property
    def time_freq_support(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._time_frequencies

    def __str__(self):
        """Get the string representation of this class."""
        txt = (
            "%s." % re.sub(r"(?<!^)(?=[A-Z])", " ", type(self).__name__)
            + "\n\n\nData Sources\n------------------------------\n"
        )
        ds_str = self._model._data_sources.__str__()
        txt += ds_str
        txt += "\n\n"
        txt += self._model.__str__()
        return txt

    def _build_op_names_from_components(self, op_name_base, components):
        # Create operator internal names based on components
        op_names = []
        if components is None:
            op_names = [op_name_base]
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if isinstance(components, str) or isinstance(components, int):
                raise ValueError("Argument 'components' must be a list.")
            for comp in components:
                if not (isinstance(comp, str) or isinstance(comp, int)):
                    raise ValueError(
                        "Argument 'components' must be a list of integers and/or strings."
                    )
                if isinstance(comp, int):
                    comp = str(comp)
                if comp not in self._component_id_to_str.keys():
                    raise ValueError(
                        f"Component {comp} is not valid. Please use one of: "
                        f"{list(self._component_id_to_str.keys())}."
                    )
                op_comp = self._component_id_to_str[comp]
                op_names.append(op_name_base + op_comp)

        # Take unique values
        return list(set(op_names))

    def _build_op_names_from_principal_components(self, op_name_base, components):
        # Create operator internal names based on principal components
        op_names = []
        if components is None:
            op_names = [op_name_base + "1", op_name_base + "2", op_name_base + "3"]
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if isinstance(components, str) or isinstance(components, int):
                raise ValueError("Argument 'components' must be a list.")
            for comp in components:
                if isinstance(comp, int):
                    comp = str(comp)
                if comp not in ["1", "2", "3"]:
                    raise ValueError("A principal component must be 1, 2, or 3.")
                op_names.append(op_name_base + comp)

        # Take unique values
        return list(set(op_names))


class MechanicalSimulation(Simulation, ABC):
    """Base class for mechanical type simulations.

    This class provides common methods and properties for all mechanical type simulations.
    """

    def __init__(self, data_sources: core.DataSources, model: core.Model):
        """Instantiate a mechanical type simulation."""
        super().__init__(data_sources, model)

    def _build_mesh_scoping(
        self,
        selection=None,
        nodes=None,
        elements=None,
        named_selection=None,
        location=core.locations.nodal,
    ):
        if (nodes is not None or elements is not None) and named_selection is not None:
            raise ValueError(
                "nodes/elements and named_selection are mutually exclusive"
            )

        if selection is not None and (
            nodes is not None or named_selection is not None or elements is not None
        ):
            raise ValueError(
                "selection and nodes/elements/named_selection are mutually exclusive"
            )

        mesh_scoping = None
        if selection:
            mesh_scoping = selection.mesh_scoping

        if named_selection:
            mesh_scoping = core.mesh_scoping_factory.named_selection_scoping(
                named_selection, server=self._model._server, model=self._model
            )

        if nodes:
            mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
                nodes, server=self._model._server
            )

        if elements:
            mesh_scoping = core.mesh_scoping_factory.elemental_scoping(
                element_ids=elements, server=self._model._server
            )

        if (
            location == core.locations.nodal
            and mesh_scoping.location != core.locations.nodal
        ) or (
            location == core.locations.elemental
            and mesh_scoping.location != core.locations.elemental
        ):
            mesh_scoping = core.operators.scoping.transpose(
                mesh_scoping=mesh_scoping,
                meshed_region=self.mesh._meshed_region,
                inclusive=1,
            ).eval()

        return mesh_scoping


class StaticMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical static simulations."""

    def _build_time_freq_scoping(
        self,
        selection: Union[Selection, None],
        times: Union[float, List[float], None],
        set_ids: Union[int, List[int], None],
        load_steps: Union[int, List[int], None],
        sub_steps: Union[int, List[int], None],
    ) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        # create from selection in priority
        if selection:
            return selection.time_freq_selection._evaluate_on(simulation=self)
        # else from set_ids
        if set_ids:
            if isinstance(set_ids, int):
                set_ids = [set_ids]
            return core.time_freq_scoping_factory.scoping_by_sets(
                cumulative_sets=set_ids, server=self._model._server
            )
        # else from times
        if times:
            if isinstance(times, float):
                times = [times]
            raise NotImplementedError
        # else from sub_steps and load_steps
        if sub_steps:
            if isinstance(sub_steps, int):
                sub_steps = [sub_steps]
            if isinstance(load_steps, int):
                load_steps = [load_steps]
            elif (
                not load_steps
                or (isinstance(load_steps, list) and len(load_steps)) != 1
            ):
                raise ValueError(
                    "Argument sub_steps requires argument load_steps to have one value."
                )
            # Translate to cumulative indices (set IDs)
            set_ids = []
            for sub_step in sub_steps:
                set_id = (
                    self._model.metadata.time_freq_support.get_cumulative_index(
                        step=load_steps[0] - 1, substep=sub_step
                    )
                    + 2
                )
                set_ids.append(set_id)
            return core.time_freq_scoping_factory.scoping_by_sets(
                cumulative_sets=set_ids, server=self._model._server
            )
        # else load_steps only
        if load_steps:
            if isinstance(load_steps, int):
                load_steps = [load_steps]
            return core.time_freq_scoping_factory.scoping_by_load_steps(
                load_steps=load_steps, server=self._model._server
            )
        # Otherwise, no argument was given, create a time_freq_scoping of the whole results
        return core.time_freq_scoping_factory.scoping_on_all_time_freqs(self._model)

    def _get_result(
        self,
        base_name: str,
        location: str,
        category: str,
        components: Union[str, List[str], int, List[int]],
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            base_name:
                Base name for the requested result.
            location:
                Location requested.
            category:
                Type of result requested. Can be "component", "principal", or "equivalent".
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        # Build the targeted time scoping
        time_scoping = self._build_time_freq_scoping(
            selection, times, set_ids, load_steps, sub_steps
        )

        # Build the targeted mesh scoping
        mesh_scoping = self._build_mesh_scoping(
            selection,
            nodes,
            elements,
            named_selection,
            location=location,
        )

        # Build the list of required operators
        if category == "components":
            op_names = self._build_op_names_from_components(
                op_name_base=base_name, components=components
            )
        elif category == "principal":
            op_names = self._build_op_names_from_principal_components(
                op_name_base=base_name, components=components
            )

        # Initialize a workflow
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False

        # If more than one operator is required for extraction,
        # a merging step using utility.merge_fields_containers is needed in the workflow.
        assemble_op = self._model.operator(name="merge::fields_container")
        # assemble_op = self._model.operator(name="utility::assemble_scalars_to_vectors_fc")
        wf.add_operator(operator=assemble_op)

        # Set the global output of the workflow
        wf.set_output_name("out", assemble_op.outputs.merged_fields_container)

        # For each required operator
        for pin, op_name in enumerate(op_names):
            # Instantiate the operator
            op = self._model.operator(name=op_name)
            # Set the time_scoping if necessary
            if time_scoping:
                op.connect(0, time_scoping)
            # Set the mesh_scoping if necessary
            if mesh_scoping:
                op.connect(1, mesh_scoping)

            op.connect(9, location)

            # Connect its output to the merge operator
            assemble_op.connect(pin=pin, inpt=op.outputs.fields_container)
            wf.add_operator(operator=op)

        return DataObject(
            wf.get_output("out", core.types.fields_container),
            columns=op_names,
            mesh_scoping=mesh_scoping,
        )

    def displacement(
        self,
        components: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract displacement results from the simulation.

        Args:
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                Times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                Load steps to get results for.
            sub_steps:
                Sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements whose nodes to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="U",
            location=core.locations.nodal,
            category="components",
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress(
        self,
        components: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental_nodal,
            category="components",
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_stress(
        self,
        components: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental,
            category="components",
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def nodal_stress(
        self,
        components: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.nodal,
            category="components",
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_principal_stress(
        self,
        components: Union[List[str], List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[List[int], None] = None,
        load_steps: Union[List[int], None] = None,
        sub_steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.elemental,
            category="principal",
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def nodal_principal_stress(
        self,
        components: Union[List[str], List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[List[int], None] = None,
        load_steps: Union[List[int], None] = None,
        sub_steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            components:
                Components to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=core.locations.nodal,
            category="principal",
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_von_mises_eqv_stress(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[List[int], None] = None,
        load_steps: Union[List[int], None] = None,
        sub_steps: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            load_steps:
                List of load steps to get results for.
            sub_steps:
                List of sub-steps to get results for. Requires load_steps to be defined.
            elements:
                List of elements to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        # Build the targeted time scoping
        time_scoping = self._build_time_freq_scoping(
            selection, times, set_ids, load_steps, sub_steps
        )

        # Build the targeted mesh scoping
        mesh_scoping = self._build_mesh_scoping(
            selection, None, elements, named_selection
        )

        # Instantiate the operator
        op = self._model.operator(name="S_eqv")
        # Set the time_scoping if necessary
        if time_scoping:
            op.connect(0, time_scoping)
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            op.connect(1, mesh_scoping)

        op.connect(9, "Elemental")

        # Initialize a workflow
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False
        wf.add_operator(operator=op)

        # Set the global output of the workflow
        wf.set_output_name("out", op.outputs.fields_container)

        return DataObject(
            wf.get_output("out", core.types.fields_container),
            columns="S_VM",
            mesh_scoping=mesh_scoping,
        )


class TransientMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical transient simulations."""


class ModalMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical modal simulations."""


class HarmonicMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical harmonic simulations."""
