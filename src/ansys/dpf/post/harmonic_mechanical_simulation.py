"""Module containing the ``HarmonicMechanicalSimulation`` class.

HarmonicMechanicalSimulation
----------------------------

"""
from typing import List, Tuple, Union
import warnings

from ansys.dpf import core as dpf
from ansys.dpf.post import locations
from ansys.dpf.post.dataframe import DataFrame
from ansys.dpf.post.index import (
    CompIndex,
    LabelIndex,
    MeshIndex,
    MultiIndex,
    ResultsIndex,
    SetIndex,
)
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
        sweeping_phase: Union[float, None] = 0.0,
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
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract results from the simulation.

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
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            category:
                Type of result requested. See the :class:`ResultCategory` class.
            components:
                Components to get results for.
            norm:
                Whether to return the norm of the results.
            amplitude:
                Whether to return the amplitude of the result. Overrides `sweeping_phase`.
            sweeping_phase:
                Sweeping phase value to extract result for. If a single `float` value is given, it
                must be in degrees. Unused if `amplitude=True`.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            node_ids:
                List of IDs of nodes to get results for.
            element_ids:
                List of IDs of elements to get results for.
            named_selections:
                Named selection or list of named selections to get results for.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

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

        comp, to_extract, columns = self._create_components(
            base_name, category, components
        )

        # Initialize a workflow
        wf = dpf.Workflow(server=self._model._server)
        wf.progress_bar = False

        if category == ResultCategory.equivalent and base_name[0] == "E":
            force_elemental_nodal = True
        # elif sweeping_phase is not None:
        #     force_elemental_nodal = True
        else:
            force_elemental_nodal = False

        # Instantiate the main result operator
        result_op = self._build_result_operator(
            name=base_name,
            location=location,
            force_elemental_nodal=force_elemental_nodal,
        )

        # Treat cyclic cases
        result_op = self._treat_cyclic(expand_cyclic, phase_angle_cyclic, result_op)

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
            base_name += "_VM"

        # Add an optional component selection step if result is vector, or matrix
        if (
            category
            in [
                ResultCategory.vector,
                ResultCategory.matrix,
            ]
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
            if len(to_extract) == 1:
                base_name += f"_{comp[0]}"
                comp = None

        # Add an optional sweeping phase or amplitude operation if requested
        # (must be after comp_selector for U)
        # (must be before norm operation for U)
        if sweeping_phase is not None and not amplitude:
            if isinstance(sweeping_phase, int):
                sweeping_phase = float(sweeping_phase)
            if not isinstance(sweeping_phase, float):
                raise ValueError("Argument sweeping_phase must be a float.")
            sweeping_op = self._model.operator(name="sweeping_phase_fc")
            sweeping_op.connect(0, out)
            sweeping_op.connect(2, sweeping_phase)
            sweeping_op.connect(3, "degree")
            sweeping_op.connect(4, False)
            wf.add_operator(operator=sweeping_op)
            out = sweeping_op.outputs.fields_container
        elif amplitude:
            amplitude_op = self._model.operator(name="amplitude_fc")
            amplitude_op.connect(0, out)
            wf.add_operator(operator=amplitude_op)
            out = amplitude_op.outputs.fields_container

        # Add an optional norm operation if requested
        # (must be after sweeping_phase for U)
        if norm:
            norm_op = self._model.operator(name="norm_fc")
            norm_op.connect(0, out)
            wf.add_operator(operator=norm_op)
            out = norm_op.outputs.fields_container
            comp = None
            base_name += "_N"

        # Set the workflow output
        wf.set_output_name("out", out)
        # Evaluate  the workflow
        fc = wf.get_output("out", dpf.types.fields_container)

        disp_wf = self._generate_disp_workflow(fc, selection)

        # Test for empty results
        if (len(fc) == 0) or all([len(f) == 0 for f in fc]):
            warnings.warn(
                message=f"Returned Dataframe with columns {columns} is empty.",
                category=UserWarning,
            )
        comp_index = None
        if comp is not None:
            comp_index = CompIndex(values=comp)
        row_indexes = [MeshIndex(location=location, fc=fc)]
        if comp_index is not None:
            row_indexes.append(comp_index)
        column_indexes = [
            ResultsIndex(values=[base_name]),
            SetIndex(values=fc.get_available_ids_for_label("time")),
        ]
        label_indexes = []
        for label in fc.labels:
            if label not in ["time"]:
                label_indexes.append(
                    LabelIndex(name=label, values=fc.get_available_ids_for_label(label))
                )

        column_indexes.extend(label_indexes)
        column_index = MultiIndex(indexes=column_indexes)

        row_index = MultiIndex(
            indexes=row_indexes,
        )

        # Return the result wrapped in a DPF_Dataframe
        df = DataFrame(
            data=fc,
            columns=column_index,
            index=row_index,
        )
        df._disp_wf = disp_wf
        return df

    def displacement(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        # amplitude: bool = False,
        # sweeping_phase: Union[float, None] = 0.0,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
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
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="U",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def velocity(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        # amplitude: bool = False,
        # sweeping_phase: Union[float, None] = 0.0,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
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
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def acceleration(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        # amplitude: bool = False,
        # sweeping_phase: Union[float, None] = 0.0,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
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
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="A",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[List[str], List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract principal stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=location,
            category=ResultCategory.principal,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[List[str], List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[List[str], List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal principal stress results from the simulation.

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
                List of IDs pf elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_eqv_von_mises(
        self,
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
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract equivalent von Mises stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_eqv_von_mises_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def stress_eqv_von_mises_nodal(
        self,
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
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal equivalent von Mises stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract stress results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are "X", "Y", "Z", "XX", "XY",
                "XZ", and their respective equivalents 1, 2, 3, 4, 5, 6.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract principal elastic strain results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=location,
            category=ResultCategory.principal,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal principal elastic strain results from the simulation.

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
                List of IDs pf elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            components:
                Components to get results for. Available components are: 1, 2, and 3.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_eqv_von_mises(
        self,
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
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract equivalent von Mises elastic strain results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_eqv_von_mises_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elastic_strain_eqv_von_mises_nodal(
        self,
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
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal equivalent von Mises elastic strain results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPEL",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def reaction_force(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract reaction force results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="RF",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elemental_volume(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental volume results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="ENG_VOL",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def elemental_mass(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental mass results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="ElementalMass",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_nodal_forces(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
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
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract element nodal forces results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=location,
            category=ResultCategory.vector,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_nodal_forces_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract element nodal forces nodal results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_nodal_forces_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract element nodal forces elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="ENF",
            location=locations.elemental,
            category=ResultCategory.vector,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def nodal_force(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            element_ids:
                List of IDs of elements to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="F",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def nodal_moment(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal moment results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            element_ids:
                List of IDs of elements to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="M",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            amplitude=False,
            sweeping_phase=None,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_centroids(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract element centroids results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="centroids",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def thickness(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract element thickness results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="thickness",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_orientations(
        self,
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
        location: Union[locations, str] = locations.elemental_nodal,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract element orientations results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=location,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_orientations_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract elemental element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            element_ids:
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

    def element_orientations_nodal(
        self,
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
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> DataFrame:
        """Extract nodal element orientations results from the simulation.

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
                List of IDs of elements to get results for.
            frequencies:
                Frequency value or list of frequency values to get results for.
            set_ids:
                Sets to get results for.
                A set is defined as a unique combination of {time, load step, sub-step}.
            all_sets:
                Whether to get results for all sets.
            load_steps:
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.
            expand_cyclic:
                For cyclic problems, whether to expand the sectors.
                Can take a list of sector numbers to select specific sectors to expand
                (one-based indexing).
                If the problem is multi-stage, can take a list of lists of sector numbers, ordered
                by stage.
            phase_angle_cyclic:
                 For cyclic problems, phase angle to apply (in degrees).

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EUL",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            amplitude=False,
            sweeping_phase=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )
