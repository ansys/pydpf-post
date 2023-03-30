"""Module containing the ``FluidSimulation`` class.

FluidSimulation
---------------

"""
from os import PathLike
from typing import List, Tuple, Union

from ansys.dpf import core as dpf
from ansys.dpf.post import locations
from ansys.dpf.post.dataframe import DataFrame
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import ResultCategory, Simulation

# Load the CFF plugin
dpf.core.load_library(
    r"D:\ANSYSDev\dpf_standalone\cff_fixes\v232\ansys\dpf\server_2023_2_pre1\dpf\plugins\Ans.Dpf.CFF.dll",  # noqa
    "cff_ops",
)


class FluidSimulation(Simulation):
    """Base class for fluid type simulations.

    This class provides common methods and properties for all fluid type simulations.
    """

    def __init__(self, result_file: Union[PathLike, str]):
        """Instantiate a mechanical type simulation."""
        model = dpf.Model(result_file)
        data_sources = model.metadata.data_sources
        super().__init__(data_sources=data_sources, model=model)

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
        # # Build the targeted time scoping
        # tot = (
        #     (set_ids is not None)
        #     + (all_sets is True)
        #     + (times is not None)
        #     + (load_steps is not None)
        #     + (selection is not None)
        # )
        # if tot > 1:
        #     raise ValueError(
        #         "Arguments all_sets, selection, set_ids, times, "
        #         "and load_steps are mutually exclusive."
        #     )

        # selection = self._build_selection(
        #     selection=selection,
        #     set_ids=set_ids,
        #     times=times,
        #     load_steps=load_steps,
        #     all_sets=all_sets,
        #     node_ids=node_ids,
        #     element_ids=element_ids,
        #     named_selections=named_selections,
        #     location=location,
        # )

        comp, to_extract, columns = self._create_components(
            base_name, category, components
        )

        # Initialize a workflow
        wf = dpf.Workflow(server=self._model._server)
        wf.progress_bar = False

        # Instantiate the main result operator
        result_op = self._build_result_operator(
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
        out = getattr(result_op.outputs, base_name.split("::")[-1])
        # # Its inputs are selected as workflow inputs for merging with selection workflows
        # wf.set_input_name("time_scoping", result_op.inputs.time_scoping)
        # wf.set_input_name("mesh_scoping", result_op.inputs.mesh_scoping)
        #
        # wf.connect_with(
        #     selection.time_freq_selection._selection,
        #     output_input_names=("scoping", "time_scoping"),
        # )
        # wf.connect_with(
        #     selection.spatial_selection._selection,
        #     output_input_names=("scoping", "mesh_scoping"),
        # )

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
            base_name="cff::cas::sv_density",
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
