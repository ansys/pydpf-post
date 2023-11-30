"""Module containing the ``ModalMechanicalSimulation`` class.

ModalMechanicalSimulation
-------------------------

"""
from typing import List, Union

from ansys.dpf import core as dpf
from ansys.dpf.post import locations
from ansys.dpf.post.dataframe import DataFrame
from ansys.dpf.post.selection import Selection, _WfNames
from ansys.dpf.post.simulation import MechanicalSimulation, ResultCategory


class ModalMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical modal simulations."""

    def _get_result_workflow(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
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
            if len(to_extract) == 1:
                base_name += f"_{comp[0]}"
                comp = None

        # Add an optional norm operation if requested
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
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        modes: Union[int, List[int], None] = None,
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
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        frequencies:
            Frequency value or list of frequency values to get results for.
        set_ids:
            List of sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets/modes.
        modes:
            Mode number or list of mode numbers to get results for.
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
            + (modes is not None)
            + (selection is not None)
        )
        if tot > 1:
            raise ValueError(
                "Arguments all_sets, selection, set_ids, frequencies, "
                "and modes are mutually exclusive."
            )
        elif tot == 0:
            set_ids = 1

        selection = self._build_selection(
            base_name=base_name,
            category=category,
            selection=selection,
            set_ids=set_ids if set_ids else modes,
            times=frequencies,
            load_steps=None,
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
            selection=selection,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
        )

        # Evaluate  the workflow
        fc = wf.get_output("out", dpf.types.fields_container)

        disp_wf = self._generate_disp_workflow(fc, selection)

        submesh = None
        if selection.outputs_mesh:
            selection.spatial_selection._selection.progress_bar = False
            submesh = selection.spatial_selection._selection.get_output(
                _WfNames.mesh, dpf.types.meshed_region
            )

        _, _, columns = self._create_components(base_name, category, components)
        return self._create_dataframe(
            fc, location, columns, comp, base_name, disp_wf, submesh
        )

    def displacement(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract displacement results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.ModalMechanicalSimulation(examples.download_modal_frame())
        >>> # Extract the displacement field for the first mode
        >>> displacement = simulation.displacement()
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (mm)
                     set_ids           1
         node_ids components
              367          X -2.9441e-01
                           Y  1.2397e+00
                           Z  5.1160e-01
              509          X -3.4043e-01
                           Y  1.8414e+00
                           Z  3.4187e-01
              ...        ...         ...
        >>> # Extract the displacement field for the first two modes
        >>> displacement = simulation.displacement(modes=[1, 2])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (mm)
                     set_ids           1           2
         node_ids components
              367          X -2.9441e-01  1.7382e+00
                           Y  1.2397e+00  5.4243e-01
                           Z  5.1160e-01 -4.2969e-01
              509          X -3.4043e-01  2.4632e+00
                           Y  1.8414e+00  7.5043e-01
                           Z  3.4187e-01 -2.7130e-01
              ...        ...         ...         ...
        >>> # Extract the displacement field for all modes
        >>> displacement = simulation.displacement(all_sets=True)
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (mm)
                     set_ids           1           2           3           4           5           6
         node_ids components
              367          X -2.9441e-01  1.7382e+00 -1.0401e-01 -3.6455e-01  2.8577e+00 -6.7501e-01
                           Y  1.2397e+00  5.4243e-01  1.9069e+00  2.1373e+00 -5.0887e-02 -1.0978e+00
                           Z  5.1160e-01 -4.2969e-01  6.5813e-01  6.7056e-01 -8.8191e-01 -1.4610e-01
              509          X -3.4043e-01  2.4632e+00 -3.1666e-01 -3.1348e-01  3.9674e+00 -5.1783e-01
                           Y  1.8414e+00  7.5043e-01  2.5367e+00  3.0538e+00 -6.2025e-02 -1.1483e+00
                           Z  3.4187e-01 -2.7130e-01  4.4146e-01  3.9606e-01 -5.0972e-01 -1.1397e-01
              ...        ...         ...         ...         ...         ...         ...         ...
        >>> # Extract the norm of the displacement field for the first mode
        >>> displacement = simulation.displacement(norm=True)
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
          results   U_N (mm)
          set_ids          1
         node_ids
              367 1.3730e+00
              509 1.9036e+00
              428 1.0166e+00
              510 1.0461e+00
             3442 1.6226e+00
             3755 1.4089e+00
              ...        ...
        >>> # Extract the displacement field along X for the first mode
        >>> displacement = simulation.displacement(components=["X"])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
          results    U_X (mm)
          set_ids           1
         node_ids
              367 -2.9441e-01
              509 -3.4043e-01
              428 -1.1434e-01
              510 -2.0561e-01
             3442 -3.1765e-01
             3755 -2.2155e-01
              ...         ...
        >>> # Extract the displacement field at nodes 23 and 24 for the first mode
        >>> displacement = simulation.displacement(node_ids=[23, 24])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (mm)
                     set_ids           1
        node_ids  components
              23           X -0.0000e+00
                           Y -0.0000e+00
                           Z -0.0000e+00
              24           X  2.8739e-02
                           Y  1.3243e-01
                           Z  1.4795e-01
        >>> # Extract the displacement field at nodes of element 40 for the first mode
        >>> displacement = simulation.displacement(element_ids=[40])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (mm)
                     set_ids           1
         node_ids components
              344          X -2.0812e-01
                           Y  1.1289e+00
                           Z  3.5111e-01
              510          X -2.0561e-01
                           Y  9.8847e-01
                           Z  2.7365e-01
              ...        ...         ...
        >>> # For cyclic results
        >>> simulation = post.ModalMechanicalSimulation(examples.find_simple_cyclic())
        >>> # Extract the displacement field with cyclic expansion on all sectors at first mode
        >>> displacement = simulation.displacement(expand_cyclic=True)
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (m)
                     set_ids          1
         node_ids components
                1          X 1.7611e-13
                           Y 8.5207e+01
                           Z 3.1717e-12
               52          X 2.3620e-12
                           Y 8.5207e+01
                           Z 2.1160e-12
              ...        ...        ...
        >>> # Extract the displacement field without cyclic expansion at first mode
        >>> displacement = simulation.displacement(expand_cyclic=False)
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                      results       U (m)
                      set_ids           1
                  base_sector           1
         node_ids  components
                1           X  4.9812e-13
                            Y  2.4100e+02
                            Z  8.9709e-12
               14           X -1.9511e-12
                            Y  1.9261e+02
                            Z  5.0359e-12
              ...         ...         ...
        >>> # Extract the displacement field with cyclic expansion on selected sectors at first mode
        >>> displacement = simulation.displacement(expand_cyclic=[1, 2, 3])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results      U (m)
                     set_ids          1
         node_ids components
                1          X 1.7611e-13
                           Y 8.5207e+01
                           Z 3.1717e-12
               52          X 2.3620e-12
                           Y 8.5207e+01
                           Z 2.1160e-12
              ...        ...        ...
        >>> # For multi-stage cyclic results
        >>> simulation = post.ModalMechanicalSimulation(
        ...     examples.download_multi_stage_cyclic_result()
        ... )
        >>> # Extract the displacement field with cyclic expansion on the first four sectors of the
        >>> # first stage at first mode
        >>> displacement = simulation.displacement(expand_cyclic=[1, 2, 3, 4])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results       U (m)
                     set_ids           1
         node_ids components
             1376          X  4.3586e-02
                           Y -3.0071e-02
                           Z -9.4850e-05
             4971          X  4.7836e-02
                           Y  2.2711e-02
                           Z -9.4850e-05
              ...        ...         ...
        >>> # Extract the displacement field with cyclic expansion on the first two sectors of both
        >>> # stages at first mode
        >>> displacement = simulation.displacement(expand_cyclic=[[1, 2], [1, 2]])
        >>> print(displacement)  # doctest: +NORMALIZE_WHITESPACE
                     results       U (m)
                     set_ids           1
         node_ids components
             1376          X  4.3586e-02
                           Y -3.0071e-02
                           Z -9.4850e-05
             4971          X  4.7836e-02
                           Y  2.2711e-02
                           Z -9.4850e-05
              ...        ...         ...
        """
        return self._get_result(
            base_name="U",
            location=locations.nodal,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal principal stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal equivalent von Mises stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental principal elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal equivalent von Mises elastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_state_variable(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract plastic state variable results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="ENL_PSV",
            location=location,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_state_variable_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental plastic state variable results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="ENL_PSV",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_state_variable_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal plastic state variable results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="ENL_PSV",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=location,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=locations.nodal,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=locations.elemental,
            category=ResultCategory.matrix,
            components=components,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_principal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract principal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=location,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_principal_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal principal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=locations.nodal,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_principal_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental principal plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=locations.elemental,
            category=ResultCategory.principal,
            components=components,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=None,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_eqv(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract equivalent plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=location,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_eqv_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal equivalent plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=locations.nodal,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def plastic_strain_eqv_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental equivalent plastic strain results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="EPPL",
            location=locations.elemental,
            category=ResultCategory.equivalent,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=None,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract reaction force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental volume results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental mass results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=None,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element centroids results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element thickness results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental nodal element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract elemental element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal element orientations results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def hydrostatic_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract hydrostatic pressure element nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="ENL_HPRES",
            location=location,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def hydrostatic_pressure_nodal(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract hydrostatic pressure nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="ENL_HPRES",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )

    def hydrostatic_pressure_elemental(
        self,
        element_ids: Union[List[int], None] = None,
        frequencies: Union[float, List[float], None] = None,
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract hydrostatic pressure elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

        Arguments `selection`, `named_selections`, and `element_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Parameters
        ----------
        element_ids:
            List of IDs of elements to get results for.
        frequencies:
            Frequency value or list of frequency values to get results for.
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            base_name="ENL_HPRES",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        location: Union[locations, str] = locations.elemental_nodal,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element nodal forces results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
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
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element nodal forces nodal results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract element nodal forces elemental results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal force results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
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
        modes: Union[int, List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
        phase_angle_cyclic: Union[float, None] = None,
        external_layer: Union[bool, List[int]] = False,
        skin: Union[bool, List[int]] = False,
    ) -> DataFrame:
        """Extract nodal moment results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `frequencies`, and `modes` are mutually
        exclusive.
        If none of the above is given, only the first mode will be returned.

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
        modes:
            Mode number or list of mode numbers to get results for.
        named_selections:
            Named selection or list of named selections to get results for.
        selection:
            Selection to get results for.
            A Selection defines both spatial and time-like criteria for filtering.
        set_ids:
            Sets to get results for. Equivalent to modes.
            Common to all simulation types for easier scripting.
        all_sets:
            Whether to get results for all sets/modes.
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
            norm=norm,
            selection=selection,
            frequencies=frequencies,
            set_ids=set_ids,
            all_sets=all_sets,
            modes=modes,
            node_ids=node_ids,
            element_ids=element_ids,
            named_selections=named_selections,
            expand_cyclic=expand_cyclic,
            phase_angle_cyclic=phase_angle_cyclic,
            external_layer=external_layer,
            skin=skin,
        )
