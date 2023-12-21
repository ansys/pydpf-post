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
from ansys.dpf.post.selection import Selection, _WfNames
from ansys.dpf.post.simulation import MechanicalSimulation, ResultCategory


class HarmonicMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical harmonic simulations."""

    def _get_result_workflow(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        amplitude: bool = False,
        sweeping_phase: Union[float, None] = 0.0,
        selection: Union[Selection, None] = None,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
    ) -> (dpf.Workflow, Union[str, list[str], None], str):
        """Generate (without evaluating) the Workflow to extract results."""
        comp, to_extract, _ = self._create_components(base_name, category, components)

        force_elemental_nodal = self._requires_manual_averaging(
            base_name=base_name,
            location=location,
            category=category,
            selection=selection,
        )

        # Instantiate the main result operator
        wf, result_op = self._build_result_workflow(
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
        if selection.requires_mesh:
            # wf.set_input_name(_WfNames.mesh, result_op.inputs.mesh)
            mesh_wf = dpf.Workflow(server=self._model._server)
            mesh_wf.add_operator(self._model.metadata.mesh_provider)
            mesh_wf.set_output_name(
                _WfNames.initial_mesh, self._model.metadata.mesh_provider
            )
            selection.spatial_selection._selection.connect_with(
                mesh_wf,
                output_input_names={_WfNames.initial_mesh: _WfNames.initial_mesh},
            )

        wf.connect_with(
            selection.spatial_selection._selection,
            output_input_names={
                "scoping": "mesh_scoping",
            },
        )

        # Treat cyclic cases
        wf = self._treat_cyclic(expand_cyclic, phase_angle_cyclic, wf)

        # Connect data_sources and streams_container inputs of selection if necessary
        if "streams" in wf.input_names:
            wf.connect("streams", self._model.metadata.streams_provider)
        if "data_sources" in wf.input_names:
            wf.connect("data_sources", self._model.metadata.data_sources)

        average_op = None
        if force_elemental_nodal:
            average_op = self._create_averaging_operator(
                location=location, selection=selection
            )

        # Add a step to compute principal invariants if result is principal
        if category == ResultCategory.principal:
            # Instantiate the required operator
            principal_op = self._model.operator(name="invariants_fc")
            # Corresponds to scripting name principal_invariants
            if average_op is not None:
                average_op[0].connect(0, out)
                principal_op.connect(0, average_op[1])
                wf.add_operators(list(average_op))
                # Set as future output of the workflow
                average_op = None
            else:
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
            wf.add_operator(operator=equivalent_op)
            # If a strain result, change the location now
            if (
                average_op is not None
                and category == ResultCategory.equivalent
                and base_name[0] == "E"
            ):
                equivalent_op.connect(0, out)
                average_op[0].connect(0, equivalent_op)
                wf.add_operators(list(average_op))
                # Set as future output of the workflow
                out = average_op[1].outputs.fields_container
            elif average_op is not None:
                average_op[0].connect(0, out)
                equivalent_op.connect(0, average_op[1])
                wf.add_operators(list(average_op))
                # Set as future output of the workflow
                out = equivalent_op.outputs.fields_container
            else:
                equivalent_op.connect(0, out)
                out = equivalent_op.outputs.fields_container
            average_op = None
            base_name += "_VM"

        if average_op is not None:
            average_op[0].connect(0, out)
            wf.add_operators(list(average_op))
            out = average_op[1].outputs.fields_container

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
            wf, out, comp, base_name = self._append_norm(wf, out, base_name)

        # Set the workflow output
        wf.set_output_name("out", out)
        wf.progress_bar = False

        return wf, comp, base_name

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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            base_name=base_name,
            category=category,
            selection=selection,
            set_ids=set_ids,
            times=frequencies,
            load_steps=load_steps,
            all_sets=all_sets,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            location=location,
            external_layer=external_layer,
            skin=skin,
        )

        wf, comp, base_name = self._get_result_workflow(
            base_name=base_name,
            location=location,
            category=category,
            components=components,
            norm=norm,
            amplitude=amplitude,
            sweeping_phase=sweeping_phase,
            selection=selection,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

        # Evaluate  the workflow
        fc = wf.get_output("out", dpf.types.fields_container)

        disp_wf = self._generate_disp_workflow(fc, selection)

        _, _, columns = self._create_components(base_name, category, components)

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

        if selection.outputs_mesh:
            selection.spatial_selection._selection.progress_bar = False
            submesh = selection.spatial_selection._selection.get_output(
                _WfNames.mesh, dpf.types.meshed_region
            )
            for i_field in range(len(fc)):
                bind_support_op = dpf.operators.utility.bind_support(
                    fc[i_field],
                    submesh,
                    server=fc._server,
                )
                fc.add_field(fc.get_label_space(i_field), bind_support_op.eval())

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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract displacement results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract acceleration results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract reaction force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental volume results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental mass results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element nodal forces results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element nodal forces nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element nodal forces elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal moment results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element centroids results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element thickness results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
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
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, the last frequency will be returned.

        Arguments `selection`, `named_selections`, `element_ids`, and `node_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
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
        external_layer:
             Select the external layer (last layer of solid elements under the skin)
             of the mesh for plotting and data extraction. If a list is passed, the external
             layer is computed over list of elements.
        skin:
             Select the skin (creates new 2D elements connecting the external nodes)
             of the mesh for plotting and data extraction. If a list is passed, the skin
             is computed over list of elements (not supported for cyclic symmetry). Getting the
             skin on more than one result (several time freq sets, split data...) is only
             supported starting with Ansys 2023R2.

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
            external_layer=external_layer,
            skin=skin,
        )
