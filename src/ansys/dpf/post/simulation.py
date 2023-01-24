"""Module containing the ``Simulation`` class."""
from abc import ABC, abstractmethod
import re
from typing import List, Union
import warnings

from ansys.dpf.core import DataSources, Model
from ansys.dpf.core.plotter import DpfPlotter

from ansys.dpf import core
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.selection import Selection


class Simulation(ABC):
    """Base class of all PyDPF-Post simulation types."""

    _component_to_id = {
        "1": 0,
        "2": 1,
        "3": 2,
        "4": 3,
        "5": 4,
        "6": 5,
        "X": 0,
        "Y": 1,
        "Z": 2,
        "XY": 3,
        "YZ": 4,
        "XZ": 5,
    }

    _component_names = ["X", "Y", "Z", "XX", "XY", "YZ"]
    _principal_names = ["1", "2", "3"]

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

    def _build_components_from_components(self, base_name, category, components):
        # Create operator internal names based on components
        out = []
        if components is None:
            out = None
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if not isinstance(components, list):
                raise ValueError(
                    "Argument 'components' must be an int, a str, or a list of either."
                )
            for comp in components:
                if not (isinstance(comp, str) or isinstance(comp, int)):
                    raise ValueError(
                        "Argument 'components' can only contain integers and/or strings."
                    )
                if isinstance(comp, int):
                    comp = str(comp)
                if comp not in self._component_to_id.keys():
                    raise ValueError(
                        f"Component {comp} is not valid. Please use one of: "
                        f"{list(self._component_to_id.keys())}."
                    )
                out.append(self._component_to_id[comp])

        # Take unique values and build names list
        if out is not None:
            out = list(set(out))
        if out is None and category == "vector":
            columns = [base_name + comp for comp in self._component_names[:3]]
        elif out is None and category == "matrix":
            columns = [base_name + comp for comp in self._component_names]
        else:
            columns = [base_name + self._component_names[i] for i in out]
        return out, columns

    def _build_components_from_principal(self, base_name, components):
        # Create operator internal names based on principal components
        out = []
        if components is None:
            return None
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if not isinstance(components, list):
                raise ValueError(
                    "Argument 'components' must be an int, a str, or a list of either."
                )
            for comp in components:
                if not (isinstance(comp, str) or isinstance(comp, int)):
                    raise ValueError(
                        "Argument 'components' can only contain integers and/or strings."
                    )
                if str(comp) not in self._principal_names:
                    raise ValueError(
                        "A principal component ID must be one of: "
                        f"{self._principal_names}."
                    )
                out.append(comp)

        # Take unique values
        out = list(set(out))
        # Build columns names
        columns = [base_name]
        if out is None:
            columns = [base_name + str(comp) for comp in self._principal_names]
        return out, columns

    @abstractmethod
    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass

    @abstractmethod
    def _build_mesh_scoping(self) -> core.mesh_scoping_factory.Scoping:
        """Generate a mesh_scoping from input arguments."""
        pass


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
    ) -> Union[core.mesh_scoping_factory.Scoping, None]:
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

        if mesh_scoping is None:
            return mesh_scoping

        # Transpose location if necessary
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
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
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
                Type of result requested. Can be "scalar", "vector", or "matrix".
            components:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
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

        # Build the list of requested results
        to_extract = None
        if category in ["vector", "matrix"]:
            to_extract, columns = self._build_components_from_components(
                base_name=base_name, category=category, components=components
            )
        elif category == "principal":
            to_extract, columns = self._build_components_from_principal(
                base_name=base_name, components=components
            )
        elif category == "scalar":
            columns = [base_name]
        else:
            raise ValueError(f"'{category}' is not a valid category value.")

        # Initialize a workflow
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False

        # Instantiate the result operator
        op = self._model.operator(name=base_name)
        # Set the time_scoping if necessary
        if time_scoping:
            op.connect(0, time_scoping)
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            op.connect(1, mesh_scoping)

        op.connect(7, self.mesh._meshed_region)
        op.connect(9, location)

        if to_extract is not None:
            # Extract the components
            extract_op = self._model.operator(name="component_selector_fc")
            extract_op.connect(0, op.outputs.fields_container)
            extract_op.connect(1, to_extract)
            wf.add_operator(operator=extract_op)
            # assemble_op.connect(pin=0, inpt=extract_op.outputs.fields_container)
            # Set the global output of the workflow
            out = extract_op.outputs.fields_container
        else:
            out = op.outputs.fields_container

        if norm:
            norm_op = self._model.operator(name="norm_fc")
            norm_op.connect(0, out)
            wf.add_operator(operator=norm_op)
            out = norm_op.outputs.fields_container

        wf.set_output_name("out", out)
        fc = wf.get_output("out", core.types.fields_container)

        if (len(fc) == 0) or all([len(f) == 0 for f in fc]):
            warnings.warn(
                message=f"Returned Dataframe with columns {columns} is empty.",
                category=UserWarning,
            )

        return DataObject(
            fields_container=fc,
            columns=columns,
            mesh_scoping=mesh_scoping,
        )

    def displacement(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
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
            component_ids:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
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
            category="vector",
            components=component_ids,
            norm=norm,
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
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal stress results from the simulation.

        Args:
            component_ids:
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
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental stress results from the simulation.

        Args:
            component_ids:
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
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal stress results from the simulation.

        Args:
            component_ids:
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
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_principal_elemental(
        self,
        component_ids: Union[List[str], List[int], None] = None,
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
            component_ids:
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
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_principal_nodal(
        self,
        component_ids: Union[List[str], List[int], None] = None,
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
            component_ids:
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
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_eqv_von_mises_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[List[int], None] = None,
        load_steps: Union[List[int], None] = None,
        sub_steps: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent Von Mises stress results from the simulation.

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
        return self._get_result(
            base_name="S_eqv",
            location=core.locations.elemental,
            category="scalar",
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_eqv_von_mises(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[List[int], None] = None,
        load_steps: Union[List[int], None] = None,
        sub_steps: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract elemental nodal equivalent Von Mises stress results from the simulation.

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
        return self._get_result(
            base_name="S_eqv",
            location=core.locations.elemental_nodal,
            category="scalar",
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stress_eqv_von_mises_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[List[int], None] = None,
        load_steps: Union[List[int], None] = None,
        sub_steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent Von Mises stress results from the simulation.

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
            base_name="S_eqv",
            location=core.locations.nodal,
            category="scalar",
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPEL",
            location=core.locations.elemental_nodal,
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPEL",
            location=core.locations.nodal,
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPEL",
            location=core.locations.elemental,
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_principal_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPEL",
            location=core.locations.nodal,
            category="principal",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elastic_strain_principal_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPEL",
            location=core.locations.elemental,
            category="principal",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPPL",
            location=core.locations.elemental_nodal,
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPPL",
            location=core.locations.nodal,
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def plastic_strain_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
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
            component_ids:
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
            base_name="EPPL",
            location=core.locations.elemental,
            category="matrix",
            components=component_ids,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def reaction_force(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract reaction force results from the simulation.

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
            base_name="RF",
            location=core.locations.nodal,
            category="vector",
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def elemental_volume(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract reaction force results from the simulation.

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
        return self._get_result(
            base_name="ENG_VOL",
            location=core.locations.elemental,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )

    def stiffness_matrix_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract stiffness matrix energy results from the simulation.

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
            base_name="ENG_SE",
            location=core.locations.elemental,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def artificial_hourglass_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract artificial hourglass energy results from the simulation.

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
            base_name="ENG_AHO",
            location=core.locations.elemental,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def thermal_dissipation_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract thermal dissipation energy results from the simulation.

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
            base_name="ENG_TH",
            location=core.locations.elemental,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def kinetic_energy(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract kinetic energy results from the simulation.

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
            base_name="ENG_KE",
            location=core.locations.elemental,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def structural_temperature(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract structural temperature element nodal results from the simulation.

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
            base_name="BFE",
            location=core.locations.elemental_nodal,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def structural_temperature_nodal(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract structural temperature nodal results from the simulation.

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
            base_name="BFE",
            location=core.locations.nodal,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def structural_temperature_elemental(
        self,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract structural temperature elemental results from the simulation.

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
            base_name="BFE",
            location=core.locations.elemental,
            category="scalar",
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def element_nodal_forces(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element nodal forces results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
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
            base_name="ENF",
            location=core.locations.elemental_nodal,
            category="vector",
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def element_nodal_forces_nodal(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element nodal forces nodal results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
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
            base_name="ENF",
            location=core.locations.nodal,
            category="vector",
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=nodes,
            elements=elements,
            named_selection=named_selection,
        )

    def element_nodal_forces_elemental(
        self,
        component_ids: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        load_steps: Union[int, List[int], None] = None,
        sub_steps: Union[int, List[int], None] = None,
        elements: Union[List[int], None] = None,
        named_selection: Union[str, None] = None,
    ) -> DataObject:
        """Extract element nodal forces elemental results from the simulation.

        Args:
            component_ids:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
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
            base_name="ENF",
            location=core.locations.elemental,
            category="vector",
            components=component_ids,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            load_steps=load_steps,
            sub_steps=sub_steps,
            nodes=None,
            elements=elements,
            named_selection=named_selection,
        )


class TransientMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical transient simulations."""

    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass


class ModalMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical modal simulations."""

    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass


class HarmonicMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical harmonic simulations."""

    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass
