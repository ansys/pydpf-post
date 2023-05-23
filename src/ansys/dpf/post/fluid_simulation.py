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

# Load the CFF plugin
dpf.core.load_library(
    r"D:\ANSYSDev\dpf_standalone\cff_fixes\2023.2.pre2\ansys\dpf\server_2023_2_pre2\dpf\plugins\dpf_cff\Ans.Dpf.CFF.dll",  # noqa
    "cff_ops",
)


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
        element_ids: Union[List[int], None] = None,
        node_ids: Union[List[int], None] = None,
        location: Union[locations, str] = locations.nodal,
    ) -> Selection:
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
        if selection is not None:
            return selection
        else:
            selection = Selection(server=self._model._server)

        # Create the SpatialSelection
        if named_selections:
            selection.select_named_selection(
                named_selection=named_selections, location=location
            )
        elif element_ids is not None:
            if location == locations.nodal:
                selection.select_nodes_of_elements(elements=element_ids, mesh=self.mesh)
            else:
                selection.select_elements(elements=element_ids)
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

    def __init__(self, result_file: Union[PathLike, str]):
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
        return SpeciesList()

    @property
    def phases(self):
        """Return the list of Phases in the simulation."""
        return Phases()

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
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[int], None] = None,
        named_selections: Union[List[str], str, None] = None,
    ) -> DataFrame:
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
            element_ids:
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
            element_ids=element_ids,
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
        if zone_ids:
            result_op.connect(25, zone_ids)
        if phases:
            label_space = {"phase": phases[0]}
            result_op.connect(1000, label_space)
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
        # Evaluate  the workflow
        fc = wf.get_output("out", dpf.types.fields_container)
        # print(fc)

        # disp_wf = self._generate_disp_workflow(fc, selection)

        return self._create_dataframe(
            fc, location, columns, comp, base_name.split("::")[-1], None
        )

    def density(
        self,
        node_ids: Union[List[int], None] = None,
        element_ids: Union[List[int], None] = None,
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
            element_ids=element_ids,
            zone_ids=zone_ids,
            phases=phases,
            named_selections=named_selections,
        )
