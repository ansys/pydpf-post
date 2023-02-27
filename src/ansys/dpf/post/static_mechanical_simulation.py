"""Module containing the ``StaticMechanicalSimulation`` class."""
from typing import List, Tuple, Union
import warnings

from ansys.dpf import core
from ansys.dpf.post import locations
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import MechanicalSimulation, ResultCategory


class StaticMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical static simulations."""

    def _get_result(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
    ) -> DataObject:
        """Extract results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

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
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            times:
                List of times to get results for.
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
        # Build the targeted time scoping
        tot = (
            (set_ids is not None)
            + (all_sets is True)
            + (times is not None)
            + (load_steps is not None)
            + (selection is not None)
        )
        if tot > 1:
            raise ValueError(
                "Arguments all_sets, selection, set_ids, times, "
                "and load_steps are mutually exclusive."
            )

        selection = self._build_selection(
            selection=selection,
            set_ids=set_ids,
            times=times,
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

        if category == ResultCategory.equivalent and base_name[0] == "E":
            force_elemental_nodal = True
        else:
            force_elemental_nodal = False

        # Instantiate the main result operator
        result_op = self._build_result_operator(
            name=base_name,
            location=location,
            force_elemental_nodal=force_elemental_nodal,
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
            principal_op = self._model.operator(name="invariants_fc")
            # Corresponds to scripting name principal_invariants
            principal_op.connect(0, out)
            wf.add_operator(operator=principal_op)
            # Set as future output of the workflow
            if len(to_extract) == 1:
                out = getattr(principal_op.outputs, f"fields_eig_{to_extract[0]+1}")
            else:
                raise NotImplementedError("Cannot combine principal results yet.")
                # We need to define the behavior for storing different results in a DataFrame

        # Add a step to compute equivalent if result is equivalent
        elif category == ResultCategory.equivalent:
            equivalent_op = self._model.operator(name="eqv_fc")
            equivalent_op.connect(0, out)
            wf.add_operator(operator=equivalent_op)
            # Set as future output of the workflow
            out = equivalent_op.outputs.fields_container
            # If a strain result, change the location now
            if force_elemental_nodal:
                average_op = None
                if location == locations.nodal:
                    average_op = self._model.operator(name="to_nodal_fc")
                elif location == locations.elemental:
                    average_op = self._model.operator(name="to_elemental_fc")
                if average_op is not None:
                    average_op.connect(0, out)
                    wf.add_operator(operator=average_op)
                    # Set as future output of the workflow
                    out = average_op.outputs.fields_container

        # Add an optional component selection step if result is vector, matrix, or principal
        if (category in [ResultCategory.vector, ResultCategory.matrix]) and (
            to_extract is not None
        ):
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

        extract_scoping = self._model.operator(name="extract_scoping")
        extract_scoping.connect(0, out)
        merge_scopings = self._model.operator(name="merge::scoping")
        merge_scopings.connect(0, extract_scoping.outputs.mesh_scoping_as_scoping)
        wf.set_output_name("scoping", merge_scopings.outputs.merged_scoping)

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
            index=wf.get_output("scoping", core.types.scoping).ids,
        )

    def displacement(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract displacement results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements whose nodes to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="S",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="S",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        components: Union[List[str], List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.

        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=location,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        components: Union[List[str], List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="S",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        components: Union[List[str], List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="S",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_eqv_von_mises(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal equivalent Von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="S",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_eqv_von_mises_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent Von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="S",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stress_eqv_von_mises_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent Von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="S",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPEL",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPEL",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=location,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPEL",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPEL",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_eqv_von_mises(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal equivalent Von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_eqv_von_mises_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent Von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPEL",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elastic_strain_eqv_von_mises_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent Von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPEL",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_state_variable(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal plastic state variable results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENL_PSV",
            location=location,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_state_variable_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental plastic state variable results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENL_PSV",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_state_variable_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal plastic state variable results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENL_PSV",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPPL",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPPL",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal principal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=location,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal principal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPPL",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental principal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPPL",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_eqv(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal equivalent plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EPPL",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_eqv_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPPL",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def plastic_strain_eqv_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EPPL",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ECR",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ECR",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ECR",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal principal creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ECR",
            location=location,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal principal creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ECR",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental principal creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ECR",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_eqv(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal equivalent creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ECR",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_equivalent_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal equivalent creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ECR",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def creep_strain_equivalent_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental equivalent creep strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ECR",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def reaction_force(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract reaction force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="RF",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elemental_volume(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental volume results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENG_VOL",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elemental_mass(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental mass results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ElementalMass",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def elemental_heat_generation(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental heat generation results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EHC",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_centroids(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract element centroids results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="centroids",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def thickness(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract element thickness results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="thickness",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_orientations(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract elemental nodal element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=location,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_orientations_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract elemental element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EUL",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_orientations_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="EUL",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def stiffness_matrix_energy(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract stiffness matrix energy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENG_SE",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def artificial_hourglass_energy(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract artificial hourglass energy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENG_AHO",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def thermal_dissipation_energy(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract thermal dissipation energy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENG_TH",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def kinetic_energy(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract kinetic energy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENG_KE",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def hydrostatic_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract hydrostatic pressure element nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENL_HPRES",
            location=location,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def hydrostatic_pressure_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract hydrostatic pressure nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENL_HPRES",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def hydrostatic_pressure_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract hydrostatic pressure elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENL_HPRES",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def structural_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract structural temperature element nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="BFE",
            location=location,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def structural_temperature_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract structural temperature nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="BFE",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def structural_temperature_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract structural temperature elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of times to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="BFE",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components="",
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_nodal_forces(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
    ) -> DataObject:
        """Extract element nodal forces results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in `locations`
                and are: "Nodal", "Elemental", and "ElementalNodal". The default "ElementalNodal"
                gives results with a value for every node at each element. "Elemental" gives results
                with one value for each element. "Nodal" gives results with one value for each node.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=location,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_nodal_forces_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract element nodal forces nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENF",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def element_nodal_forces_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract element nodal forces elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="ENF",
            location=locations.elemental,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def nodal_force(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements whose nodes to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="F",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )

    def nodal_moment(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataObject:
        """Extract nodal moment results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements whose nodes to get results for.
            times:
                List of time values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z",
                and their respective equivalents 1, 2, 3.
            norm:
                Whether to return the norm of the results.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load steps to get results for.
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
            base_name="M",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
        )
