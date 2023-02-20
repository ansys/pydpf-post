"""Module containing the ``HarmonicMechanicalSimulation`` class."""
from typing import List, Tuple, Union
import warnings

from ansys.dpf import core
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import MechanicalSimulation, ResultCategory


class HarmonicMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical harmonic simulations."""

    def _get_result(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        amplitude: bool = False,
        sweeping_phase: Union[float, None] = None,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

         Args:
            base_name:
                Base name for the requested result.
            location:
                Location requested.
            category:
                Type of result requested. See the :class:`ResultCategory` class.
            components:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
            amplitude:
                Whether to return the amplitude of the result.
            sweeping_phase:
                Sweeping phase value to extract result for. If a single `float` value is given, it
                must be in degrees. Mutually exclusive with `amplitude`.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                List of sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                List of load steps to get results for.
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            named_selections:
                Named selection or list of named selections to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        # Build the targeted spatial and time scoping
        tot = (
            (set_ids is not None)
            + (all_sets is True)
            + (frequencies is not None)
            + (load_steps is not None)
            + (selection is not None)
        )
        if tot > 1:
            raise ValueError(
                "Arguments all_sets, selection, set_ids, frequencies, "
                "and load_steps are mutually exclusive."
            )
        elif tot == 0:
            set_ids = 1

        tot = (
            (node_ids is not None)
            + (element_ids is not None)
            + (named_selections is not None)
            + (selection is not None)
        )
        if tot > 1:
            raise ValueError(
                "Arguments selection, named_selections, element_ids, "
                "and node_ids are mutually exclusive"
            )

        selection = self._build_selection(
            selection=selection,
            set_ids=set_ids,
            times=frequencies,
            load_steps=load_steps,
            all_sets=all_sets,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            location=location,
        )

        # Build the list of requested results
        if category in [ResultCategory.scalar, ResultCategory.equivalent]:
            # A scalar or equivalent result has no components
            to_extract = None
            columns = [base_name]
        elif category in [ResultCategory.vector, ResultCategory.matrix]:
            # A matrix or vector result can have components selected
            to_extract, columns = self._build_components_from_components(
                base_name=base_name, category=category, components=components
            )
        elif category == ResultCategory.principal:
            # A principal type of result can have components selected
            to_extract, columns = self._build_components_from_principal(
                base_name=base_name, components=components
            )
        else:
            raise ValueError(f"'{category}' is not a valid category value.")

        # Initialize a workflow
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False

        # Instantiate the main result operator
        result_op = self._build_result_operator(
            name=base_name,
            location=location,
        )
        # Its output is selected as future workflow output for now
        out = result_op.outputs.fields_container
        # Its inputs are selected as workflow inputs for merging with selection workflows
        wf.set_input_name("time_scoping", result_op.inputs.time_scoping)
        wf.set_input_name("mesh_scoping", result_op.inputs.mesh_scoping)

        wf.connect_with(
            selection.time_freq_selection._selection,
            output_input_names=("scoping", "time_scoping"),
        )
        wf.connect_with(
            selection.spatial_selection._selection,
            output_input_names=("scoping", "mesh_scoping"),
        )

        # Connect data_sources and streams_container inputs of selection if necessary
        if "streams" in wf.input_names:
            wf.connect("streams", self._model.metadata.streams_provider)
        if "data_sources" in wf.input_names:
            wf.connect("data_sources", self._model.metadata.data_sources)

        # Add a step to compute principal invariants if result is principal
        if category == ResultCategory.principal:
            # Instantiate the required operator
            principal_op = self._model.operator(name="eig_values_fc")
            principal_op.connect(0, out)
            wf.add_operator(operator=principal_op)
            # Set as future output of the workflow
            out = principal_op.outputs.fields_container

        # Add a step to compute equivalent if result is equivalent
        elif category == ResultCategory.equivalent:
            # If a stress result, use one operator
            if base_name[0] == "S":
                equivalent_op = self._model.operator(name="segalmaneqv_fc")
            # If a strain result, use another
            elif base_name[0] == "E":
                equivalent_op = self._model.operator(name="eqv_fc")
            # Throw otherwise
            else:
                raise ValueError(
                    f"Category {ResultCategory.equivalent} "
                    "is only available for stress or strain results."
                )
            equivalent_op.connect(0, out)
            wf.add_operator(operator=equivalent_op)
            # Set as future output of the workflow
            out = equivalent_op.outputs.fields_container

        # Add an optional component selection step if result is vector, matrix, or principal
        if (
            category
            in [ResultCategory.vector, ResultCategory.matrix, ResultCategory.principal]
        ) and (to_extract is not None):
            # Instantiate a component selector operator
            extract_op = self._model.operator(name="component_selector_fc")
            # Feed it the current workflow output
            extract_op.connect(0, out)
            # Feed it the requested components
            extract_op.connect(1, to_extract)
            wf.add_operator(operator=extract_op)
            # Set as future output of the workflow
            out = extract_op.outputs.fields_container

        # Add an optional norm operation if requested
        if norm:
            norm_op = self._model.operator(name="norm_fc")
            norm_op.connect(0, out)
            wf.add_operator(operator=norm_op)
            out = norm_op.outputs.fields_container

        # Set the workflow output
        wf.set_output_name("out", out)
        # Evaluate  the workflow
        fc = wf.get_output("out", core.types.fields_container)

        # Test for empty results
        if (len(fc) == 0) or all([len(f) == 0 for f in fc]):
            warnings.warn(
                message=f"Returned Dataframe with columns {columns} is empty.",
                category=UserWarning,
            )
        # Return the result wrapped in a DPF_Dataframe
        return DataObject(
            fields_container=fc,
            columns=columns,
            mesh_scoping=None,
        )

    def displacement(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        amplitude: bool = False,
        sweeping_phase: Union[float, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract displacement results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

         Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements whose nodes to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            amplitude:
                Whether to return the amplitude of the result.
            sweeping_phase:
                Sweeping phase value to extract result for. If a single `float` value is given, it
                must be in degrees. Mutually exclusive with `amplitude`.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                List of load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="U",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            amplitude=amplitude,
            sweeping_phase=sweeping_phase,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def velocity(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        amplitude: bool = False,
        sweeping_phase: Union[float, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

         Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements whose nodes to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            amplitude:
                Whether to return the amplitude of the result.
            sweeping_phase:
                Sweeping phase value to extract result for. If a single `float` value is given, it
                must be in degrees. Mutually exclusive with `amplitude`.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                List of load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="V",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            amplitude=amplitude,
            sweeping_phase=sweeping_phase,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def acceleration(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        amplitude: bool = False,
        sweeping_phase: Union[float, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract acceleration results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

         Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements whose nodes to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            amplitude:
                Whether to return the amplitude of the result.
            sweeping_phase:
                Sweeping phase value to extract result for. If a single `float` value is given, it
                must be in degrees. Mutually exclusive with `amplitude`.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                List of load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="A",
            location=core.locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            amplitude=amplitude,
            sweeping_phase=sweeping_phase,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )
