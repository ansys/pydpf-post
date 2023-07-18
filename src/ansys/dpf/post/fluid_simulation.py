"""Module containing the ``FluidSimulation`` class.

FluidSimulation
---------------

"""
from os import PathLike
from typing import List, Union

from ansys.dpf import core as dpf
from ansys.dpf.post import locations
from ansys.dpf.post.dataframe import DataFrame
from ansys.dpf.post.phase import PhasesDict
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import ResultCategory, Simulation
from ansys.dpf.post.species import SpeciesDict
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
                    raise ValueError(f"Could not find time={t} in the simulation.")
            selection.select_time_freq_sets(
                time_freq_sets=available_times_to_extract_set_ids
            )

        else:
            # Otherwise, no argument was given, create a time_freq_scoping of the last set only
            selection.select_time_freq_sets(
                time_freq_sets=[self.time_freq_support.n_sets]
            )
        return selection

    def __init__(
        self,
        result_file: Union[PathLike, str, dpf.DataSources, None] = None,
        cas: Union[PathLike, str, List[Union[PathLike, str]], None] = None,
        dat: Union[PathLike, str, List[Union[PathLike, str]], None] = None,
        flprj: Union[PathLike, str, None] = None,
        server: Union[dpf.server_types.BaseServer, None] = None,
    ):
        """Instantiate a mechanical type simulation."""
        tot = (
            (result_file is not None)
            + (cas is not None and dat is not None)
            + (flprj is not None)
        )
        if tot > 1:
            raise ValueError(
                "Argument result_file, cas and dat, and flprj are mutually exclusive."
            )
        elif tot < 1:
            raise ValueError(
                "One of result_file, cas and dat, or flprj argument must be set."
            )
        if result_file:
            ds = result_file
        else:
            ds = dpf.DataSources()
            if flprj:
                raise NotImplementedError("flprj input not accepted yet")
            if cas:
                if not isinstance(cas, list):
                    cas = [cas]
                for c in cas:
                    ds.set_result_file_path(c, "cas")
            if dat:
                if not isinstance(dat, list):
                    dat = [dat]
                for d in dat:
                    ds.add_file_path(d, "dat")
        model = dpf.Model(ds, server=server)
        data_sources = model.metadata.data_sources
        super().__init__(data_sources=data_sources, model=model)

    @property
    def zones(self) -> Zones:
        """Return the list of Zones in the simulation."""
        return Zones()

    @property
    def species(self) -> SpeciesDict:
        """Return the list of Species in the simulation."""
        return SpeciesDict(self)

    @property
    def phases(self) -> PhasesDict:
        """Return the list of PhasesDict in the simulation."""
        return PhasesDict(self)

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
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        named_selections: Union[List[str], str, None] = None,
    ) -> DataFrame:
        """Extract results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
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
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of elements to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            + (selection is not None)
        )
        if tot > 1:
            raise ValueError(
                "Arguments all_sets, selection, set_ids, and times "
                "are mutually exclusive."
            )

        selection = self._build_selection(
            base_name=base_name,
            category=category,
            selection=selection,
            set_ids=set_ids,
            times=times,
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
        if qualifiers:
            labels = list(qualifiers.keys())
            lists_labels.extend(labels)
            lists.extend([qualifiers[key] for key in labels])
        else:
            if set_ids:
                lists.append(set_ids)
                lists_labels.append("time")
            if zone_ids:
                lists.append(zone_ids)
                lists_labels.append("zone")
            if phases:
                phase_ids = []
                available_phases = self.phases
                for phase in phases:
                    phase_ids.append(available_phases[phase].id)
                lists.append(phase_ids)
                lists_labels.append("phase")
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
        # print(fc)
        if location is None and len(fc) > 0:
            location = fc[0].location
        if location == locations.elemental:
            location = "cells"
        return self._create_dataframe(
            fc, location, columns, comp, base_name.split("::")[-1], None
        )

    def density(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def density_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def density_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `face_ids`, and `cell_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def density_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract density results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def dynamic_viscosity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def dynamic_viscosity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def dynamic_viscosity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def dynamic_viscosity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract dynamic viscosity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of cell zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def enthalpy(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def enthalpy_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def enthalpy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def enthalpy_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract enthalpy results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of cell zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def entropy(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def entropy_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def entropy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def entropy_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract entropy results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def epsilon(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def epsilon_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def epsilon_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def epsilon_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract epsilon results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mach_number(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mach_number_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mach_number_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mach_number_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mach number results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_flow_rate(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_flow_rate_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_flow_rate_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_flow_rate_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass flow rate results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_fraction(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_fraction_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_fraction_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mass_fraction_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mass fraction results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_static_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_static_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_static_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean static pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            location=locations.faces,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def mean_velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract mean velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def omega(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def omega_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def omega_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def omega_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract omega results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_static_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_static_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_static_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS static pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            location=locations.faces,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def rms_velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract RMS velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def specific_heat(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def specific_heat_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def specific_heat_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def specific_heat_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract specific heat results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def static_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def static_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def static_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract static pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def superficial_velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract superficial velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def superficial_velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract superficial velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def superficial_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract superficial velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            location=locations.faces,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def superficial_velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract superficial velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def surface_heat_rate(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract surface heat rate results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def surface_heat_rate_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract surface heat rate results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def surface_heat_rate_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract surface heat rate results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def thermal_conductivity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def thermal_conductivity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def thermal_conductivity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def thermal_conductivity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract thermal conductivity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_pressure(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_pressure_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_pressure_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total pressure results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_temperature(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_temperature_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def total_temperature_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract total temperature results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_kinetic_energy_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent kinetic energy results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_viscosity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_viscosity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_viscosity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def turbulent_viscosity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract turbulent viscosity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def velocity(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract velocity results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def velocity_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract velocity results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract velocity results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            location=locations.faces,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def velocity_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract velocity results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def volume_fraction(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def volume_fraction_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def volume_fraction_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def volume_fraction_on_cells(
        self,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract volume fraction results on cells from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `cell_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=None,
            face_ids=None,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def wall_shear_stress(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract wall shear stress results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def wall_shear_stress_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract wall shear stress results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def wall_shear_stress_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract wall shear stress results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
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
            location=locations.faces,
            category=ResultCategory.vector,
            components=components,
            norm=norm,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def y_plus(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract y+ results from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def y_plus_on_nodes(
        self,
        node_ids: Union[List[int], None] = None,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract y+ results on nodes from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        node_ids:
            List of IDs of nodes to get results for.
        face_ids:
            List of IDs of faces which nodes to get results for.
        cell_ids:
            List of IDs of cells which nodes to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            node_ids=node_ids,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )

    def y_plus_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        cell_ids: Union[List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
        times: Union[float, List[float], None] = None,
        set_ids: Union[int, List[int], None] = None,
        all_sets: bool = False,
        named_selections: Union[List[str], str, None] = None,
        selection: Union[Selection, None] = None,
    ) -> DataFrame:
        """Extract y+ results on faces from the simulation.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `cell_ids`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
        cell_ids:
            List of IDs of cells which faces to get results for.
        zone_ids:
            List of IDs of zones to get results for.
        phases:
            List of IDs of phases to get results for.
        species:
            List of IDs of species to get results for.
        qualifiers:
            Dictionary of qualifier labels with associated values to get results for.
            Overrides any other qualifier argument such as `phases`, `species` or `zone_ids`.
        times:
            List of time values to get results for.
        set_ids:
            Sets to get results for.
            A set is defined as a unique combination of {time, load step, sub-step}.
        all_sets:
            Whether to get results for all sets.
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
            location=locations.faces,
            category=ResultCategory.scalar,
            components=None,
            norm=False,
            selection=selection,
            times=times,
            set_ids=set_ids,
            all_sets=all_sets,
            node_ids=None,
            face_ids=face_ids,
            cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
        )
