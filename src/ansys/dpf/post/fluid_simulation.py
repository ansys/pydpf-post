"""Module containing the ``FluidSimulation`` class.

FluidSimulation
---------------

"""
from os import PathLike
from typing import List, Tuple, Union

from ansys.dpf import core as dpf
from ansys.dpf.post import locations
from ansys.dpf.post.dataframe import DataFrame
from ansys.dpf.post.phase import Phases
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import ResultCategory, Simulation
from ansys.dpf.post.species import SpeciesList
from ansys.dpf.post.zone import Zones


class FluidSimulation(Simulation):
    """Base class for fluid type simulations.

    This class provides common methods and properties for all fluid type simulations.
    """

    def _build_selection(
        self,
        base_name: str,
        category: ResultCategory,
        selection: Union[Selection, None],
        set_ids: Union[int, List[int], None],
        times: Union[float, List[float], None],
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        cell_ids: Union[List[int], None] = None,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        location: Union[locations, str] = locations.nodal,
    ) -> Selection:
        tot = (
            (node_ids is not None)
            + (face_ids is not None)
            + (cell_ids is not None)
            + (named_selections is not None)
            + (selection is not None)
        )
        if tot > 1:
            raise ValueError(
                "Arguments selection, named_selections, cell_ids, "
                "and node_ids are mutually exclusive"
            )
        if selection is not None:
            return selection
        else:
            selection = Selection(server=self._model._server)

        # Create the SpatialSelection
        if named_selections:
            selection.select_named_selection(
                named_selection=named_selections, location=location
            )
        elif cell_ids is not None:
            if location == locations.nodal:
                selection.select_nodes_of_elements(elements=cell_ids, mesh=self.mesh)
            else:
                selection.select_elements(elements=cell_ids)
        elif face_ids is not None:
            if location == locations.nodal:
                selection.select_nodes_of_elements(elements=face_ids, mesh=self.mesh)
            else:
                selection.select_elements(elements=face_ids)
        elif node_ids is not None:
            if location != locations.nodal:
                raise ValueError(
                    "Argument 'node_ids' can only be used if 'location' "
                    "is equal to 'post.locations.nodal'."
                )
            selection.select_nodes(nodes=node_ids)

        # Create the TimeFreqSelection
        if all_sets:
            selection.time_freq_selection.select_with_scoping(
                dpf.time_freq_scoping_factory.scoping_on_all_time_freqs(self._model)
            )
        elif set_ids is not None:
            selection.select_time_freq_sets(time_freq_sets=set_ids)

        elif times is not None:
            # Check input
            if isinstance(times, list):
                if any([not isinstance(t, (float, int)) for t in times]):
                    raise ValueError("Argument times must contain numeric values only.")
            elif isinstance(times, float) or isinstance(times, int):
                times = [times]
            else:
                raise TypeError("Argument times must be a number or a list of numbers.")

            # Get the set_ids for available time values matching the requested time values.
            available_times = self.time_freq_support.time_frequencies.data
            precision = self._get_time_freq_precision
            available_times_to_extract_set_ids = []
            last_extracted_index = -1
            len_available = len(available_times)
            for t in times:
                found = False
                i = last_extracted_index + 1
                while i < len_available:
                    if abs(float(t) - available_times[i]) < precision:
                        last_extracted_index = i
                        available_times_to_extract_set_ids.append(i + 1)
                        found = True
                    i += 1
                if not found:
                    raise ValueError(
                        f"Could not find time={t}{self.units['time/frequency']} "
                        f"in the simulation."
                    )
            selection.select_time_freq_sets(
                time_freq_sets=available_times_to_extract_set_ids
            )

        # else from load_steps
        elif load_steps is not None:
            # If load_steps and sub_steps
            if len(load_steps) == 2:
                # Translate to cumulative indices (set IDs)
                set_ids = []
                sub_steps = load_steps[1]
                if not isinstance(sub_steps, list):
                    sub_steps = [sub_steps]
                set_id_0 = self._model.metadata.time_freq_support.get_cumulative_index(
                    step=load_steps[0] - 1, substep=sub_steps[0] - 1
                )
                if set_id_0 == -1:
                    raise ValueError(
                        f"Sub-step {sub_steps[0]} of load-step {load_steps[0]} "
                        f"does not exist."
                    )
                else:
                    set_id_0 += 1
                set_ids.extend([set_id_0 + i for i in range(len(sub_steps))])
                selection.select_time_freq_sets(time_freq_sets=set_ids)

            else:
                if isinstance(load_steps, int):
                    load_steps = [load_steps]
                selection.time_freq_selection.select_load_steps(load_steps=load_steps)
            return selection

        else:
            # Otherwise, no argument was given, create a time_freq_scoping of the last set only
            selection.select_time_freq_sets(
                time_freq_sets=[self.time_freq_support.n_sets]
            )
        return selection

    def __init__(self, result_file: Union[PathLike, str, dpf.DataSources]):
        """Instantiate a mechanical type simulation."""
        model = dpf.Model(result_file)
        data_sources = model.metadata.data_sources
        super().__init__(data_sources=data_sources, model=model)

    @property
    def zones(self):
        """Return the list of Zones in the simulation."""
        return Zones()

    @property
    def species(self):
        """Return the list of Species in the simulation."""
        return SpeciesList(self)

    @property
    def phases(self):
        """Return the list of Phases in the simulation."""
        return Phases(self)

    def _get_result(
        self,
        base_name: str,
        category: ResultCategory,
        location: Union[str, None] = None,
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
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
    ) -> DataFrame:
        """Extract results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            base_name:
                Base name for the requested result.
            category:
                Type of result requested. See the :class:`ResultCategory` class.
            location:
                Location to extract results at. Available locations are listed in
                class:`post.locations` and are: `post.locations.nodal`,
                `post.locations.elemental`, and `post.locations.elemental_nodal`.
                Using the default `post.locations.elemental_nodal` results in a value
                for every node at each element. Similarly, using `post.locations.elemental`
                gives results with one value for each element, while using `post.locations.nodal`
                gives results with one value for each node.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            node_ids:
                List of IDs of nodes to get results for.
            cell_ids:
                List of IDs of elements to get results for.
            named_selections:
                Named selection or list of named selections to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

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
            base_name=base_name,
            category=category,
            selection=selection,
            set_ids=set_ids,
            times=times,
            load_steps=load_steps,
            all_sets=all_sets,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            named_selections=named_selections,
            location=location,
        )

        comp, to_extract, columns = self._create_components(
            base_name, category, components
        )

        # Initialize a workflow
        wf, result_op = self._build_result_workflow(
            name=base_name,
            location=location,
            force_elemental_nodal=False,
        )
        lists = []
        lists_labels = []
        # if zone_ids is None:
        #     zone_ids = [zone.id for zone in self.zones]
        # else:
        #     result_op.connect(25, zone_ids)
        if zone_ids:
            lists.append(zone_ids)
            lists_labels.append("zone")
        # if phases is None:
        #     phases = [phase.id for phase in self.phases]
        # else:
        #     label_space = {"phase": phases[0]}
        #     result_op.connect(1000, label_space)
        if phases:
            lists.append(phases)
            lists_labels.append("phase")
        # if species is None:
        #     species = [species_i.id for species_i in self.species]
        if species:
            lists.append(species)
            lists_labels.append("species")

        if lists:
            import itertools

            for i, c in enumerate(itertools.product(*lists)):
                label_space = {}
                for j, label in enumerate(lists_labels):
                    label_space[label] = c[j]
                result_op.connect(1000 + i, label_space)

        # Its output is selected as future workflow output for now
        # print(result_op)

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

        # Add an optional norm operation if requested
        if norm:
            norm_op = self._model.operator(name="norm_fc")
            norm_op.connect(0, out)
            wf.add_operator(operator=norm_op)
            out = norm_op.outputs.fields_container

        # Set the workflow output
        wf.set_output_name("out", out)
        wf.progress_bar = False
        # Evaluate  the workflow
        fc = wf.get_output("out", dpf.types.fields_container)
        print(fc)
        if location is None and len(fc) > 0:
            location = fc[0].location
        return self._create_dataframe(
            fc, location, columns, comp, base_name.split("::")[-1], None
        )

    def density(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="RHO",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def density_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones which nodes to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="RHO",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def density_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `face_ids`, and `cell_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of face zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="RHO",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def density_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of cell zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="RHO",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def dynamic_viscosity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MU",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def dynamic_viscosity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MU",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def dynamic_viscosity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MU",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def dynamic_viscosity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of cell zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MU",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def enthalpy(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="H_S",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def enthalpy_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="H_S",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def enthalpy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="H_S",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def enthalpy_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of cell zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="H_S",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def entropy(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S_S",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def entropy_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S_S",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def entropy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S_S",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def entropy_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="S_S",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def epsilon(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPS",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def epsilon_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPS",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def epsilon_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPS",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def epsilon_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="EPS",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mach_number(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MACH",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mach_number_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MACH",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mach_number_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MACH",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mach_number_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MACH",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_flow_rate(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MDOT",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_flow_rate_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MDOT",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_flow_rate_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MDOT",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_flow_rate_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MDOT",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_fraction(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Y",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_fraction_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Y",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_fraction_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Y",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mass_fraction_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Y",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_static_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SA",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_static_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SA",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SA",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_static_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SA",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_A",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_A",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_A",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_A",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract mean velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_A",
            location=None,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract mean velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_A",
            location=None,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract mean velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_A",
            location="face",
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def mean_velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract mean velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_A",
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
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def omega(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="OME",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def omega_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="OME",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def omega_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="OME",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def omega_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="OME",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_static_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SRMS",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_static_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SRMS",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SRMS",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_static_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_SRMS",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_RMS",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_RMS",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_RMS",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_RMS",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract RMS velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_RMS",
            location=None,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract RMS velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_RMS",
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
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract RMS velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_RMS",
            location="face",
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def rms_velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract RMS velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_RMS",
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
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def specific_heat(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="CP",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def specific_heat_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="CP",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def specific_heat_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="CP",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def specific_heat_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="CP",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def static_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_S",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def static_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_S",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_S",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def static_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_S",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def superficial_velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract superficial velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_SUP",
            location=None,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def superficial_velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract superficial velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_SUP",
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
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def superficial_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract superficial velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_SUP",
            location="face",
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def superficial_velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract superficial velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V_SUP",
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
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def surface_heat_rate(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract surface heat rate results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Q",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def surface_heat_rate_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract surface heat rate results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Q",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def surface_heat_rate_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract surface heat rate results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="Q",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def thermal_conductivity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="KT",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def thermal_conductivity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="KT",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def thermal_conductivity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="KT",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def thermal_conductivity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="KT",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_TOT",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_TOT",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_TOT",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="P_TOT",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_TOT",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_TOT",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_TOT",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def total_temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TEMP_TOT",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="K",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="K",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="K",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="K",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_viscosity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MUT",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_viscosity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MUT",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_viscosity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MUT",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def turbulent_viscosity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="MUT",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V",
            location=None,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

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
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V",
            location="face",
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="V",
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
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def volume_fraction(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="VOF",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def volume_fraction_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="VOF",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def volume_fraction_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="VOF",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def volume_fraction_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="VOF",
            location=locations.elemental,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def wall_shear_stress(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract wall shear stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TAU",
            location=None,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def wall_shear_stress_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract wall shear stress results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TAU",
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
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def wall_shear_stress_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
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
    ) -> DataFrame:
        """Extract wall shear stress results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
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
                Load step number or list of load step numbers to get results for.
                One can specify sub-steps of a load step with a tuple of format:
                (load-step, sub-step number or list of sub-step numbers).
            named_selections:
                Named selection or list of named selections to get results for.
            selection:
                Selection to get results for.
                A Selection defines both spatial and time-like criteria for filtering.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="TAU",
            location="face",
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def y_plus(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract y+ results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="YPLUS",
            location=None,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def y_plus_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract y+ results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            node_ids:
                List of IDs of nodes to get results for.
            face_ids:
                List of IDs of faces which nodes to get results for.
            cell_ids:
                List of IDs of cells which nodes to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="YPLUS",
            location=locations.nodal,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )

    def y_plus_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        load_steps: Union[
            int, List[int], Tuple[int, Union[int, List[int]]], None
        ] = None,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract y+ results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, `times`, and `load_steps` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Args:
            face_ids:
                List of IDs of faces to get results for.
            cell_ids:
                List of IDs of cells which faces to get results for.
            zone_ids:
                List of IDs of zones to get results for.
            times:
                List of time values to get results for.
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        return self._get_result(
            base_name="YPLUS",
            location="face",
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            load_steps=load_steps,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )
