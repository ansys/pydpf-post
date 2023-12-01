"""Module containing the ``FluidSimulation`` class.

FluidSimulation
---------------

"""
from os import PathLike
from typing import List, Union

from ansys.dpf.core.server_types import BaseServer

from ansys.dpf import core as dpf
from ansys.dpf.post import locations
from ansys.dpf.post.dataframe import DataFrame
from ansys.dpf.post.mesh_info import FluidMeshInfo
from ansys.dpf.post.phase import PhasesDict
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.simulation import ResultCategory, Simulation
from ansys.dpf.post.species import SpeciesDict


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
            elif location == locations.faces:
                selection.select_faces_of_elements(elements=cell_ids, mesh=self.mesh)
            else:
                selection.select_elements(elements=cell_ids)
        elif face_ids is not None:
            if location == locations.nodal:
                selection.select_nodes_of_faces(faces=face_ids, mesh=self.mesh)
            else:
                selection.select_faces(faces=face_ids)
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
        server: Union[BaseServer, None] = None,
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
                ds.set_result_file_path(flprj, "flprj")
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
        self._mesh_info = None

    @property
    def mesh_info(self) -> FluidMeshInfo:
        """Return available mesh information."""
        if not self._mesh_info:
            self._mesh_info = FluidMeshInfo(self._model.metadata.mesh_info)
        return self._mesh_info

    @property
    def cell_zones(self) -> dict:
        """Return a dictionary of the cell zones in the simulation."""
        return self.mesh_info.cell_zones

    @property
    def face_zones(self) -> dict:
        """Return a dictionary of the face zones in the simulation.

        For CFX files, we gather face zones in COMPOSITE entities more physics-related.
        """
        return self.mesh_info.face_zones

    @property
    def species(self) -> SpeciesDict:
        """Return a dictionary-like object of species in the simulation."""
        return SpeciesDict(self)

    @property
    def phases(self) -> PhasesDict:
        """Return a dictionary-like object of phases in the simulation."""
        return PhasesDict(self)

    def _filter_zones(self, zone_ids: List[int], keep: locations):
        """Filter zone IDs to only keep zones of the given type."""
        if keep == locations.elemental:
            ref = set(self.cell_zones.keys())
        elif keep == locations.faces:
            ref = set(self.face_zones.keys())
        return [i for i in zone_ids if i in ref]

    def _get_result_workflow(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        components: Union[str, List[str], int, List[int], None] = None,
        norm: bool = False,
        selection: Union[Selection, None] = None,
        set_ids: Union[int, List[int], None] = None,
        zone_ids: Union[List[int], None] = None,
        phases: Union[List[Union[int, str]], None] = None,
        species: Union[List[int], None] = None,
        qualifiers: Union[dict, None] = None,
    ) -> (dpf.Workflow, Union[str, list[str], None], str):
        """Generate (without evaluating) the Workflow to extract results."""
        comp, to_extract, columns = self._create_components(
            base_name, category, components
        )

        # Initialize a workflow
        wf, result_op = self._build_result_workflow(
            name=base_name,
            location=location,
            force_elemental_nodal=False,
        )
        query_regions_meshes = False
        lists = []
        lists_labels = []
        if qualifiers:
            labels = list(qualifiers.keys())
            lists_labels.extend(labels)
            lists.extend([qualifiers[key] for key in labels])
            if "zone" in labels:
                query_regions_meshes = qualifiers["zone"]
        else:
            if set_ids:
                lists.append(set_ids)
                lists_labels.append("time")
            if zone_ids:
                lists.append(zone_ids)
                lists_labels.append("zone")
                query_regions_meshes = zone_ids
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

        if query_regions_meshes:
            # Results have been queried on regions,
            # A MeshesProvider is required to give meshes as input of the source operator
            meshes_provider_op = self._model.operator("meshes_provider")
            meshes_provider_op.connect(25, query_regions_meshes)
            result_op.connect(7, meshes_provider_op.outputs.meshes)
            wf.add_operator(meshes_provider_op)
        else:
            # Results have been queried on the whole mesh,
            # A MeshProvider is required to give the mesh as input of the source operator
            mesh_provider_op = self._model.operator("mesh_provider")
            result_op.connect(7, mesh_provider_op.outputs.mesh)
            wf.add_operator(mesh_provider_op)

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
            wf, out, comp, base_name = self._append_norm(wf, out, base_name)

        # Set the workflow output
        wf.set_output_name("out", out)
        wf.progress_bar = False

        return wf, comp, base_name

    def _get_result(
        self,
        base_name: str,
        category: ResultCategory,
        native_location: str,
        location: Union[locations, str, None] = None,
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
        integrated: Union[bool, None] = None,
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
        native_location:
            Native location of the result.
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.
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
        integrated:
            An integrated result cannot be requested on another location than its native_location.

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

        # Raise for integrated result queried on anything else than its native location
        if integrated:
            if (native_location == "Faces" and location != locations.faces) or (
                native_location == "Elements" and location != locations.elemental
            ):
                raise ValueError(
                    f"Cannot query a {native_location} integrated result on {location}."
                )

        # Define required averaging step
        averaging_op_name = None
        if native_location == "Nodes":
            if location != locations.nodal:
                averaging_op_name = "to_elemental_fc"
        elif native_location == "Elemental":
            if location == locations.faces:
                raise ValueError("Cannot query elemental results on faces.")
            elif location == locations.nodal:
                # averaging_op_name = "to_nodal_fc"
                pass  # nodal averaging seems to be automatic
        elif native_location == "Faces":
            if location == locations.elemental:
                raise ValueError("Cannot query faces results on elements.")
            elif location == locations.nodal:
                # averaging_op_name = "to_nodal_fc"
                pass  # nodal averaging seems to be automatic
        elif native_location == "ElementalAndFaces":
            if location == locations.nodal:
                # averaging_op_name = "to_nodal_fc"
                pass  # nodal averaging seems to be automatic
            elif location == locations.faces:
                if qualifiers and ("zone" in qualifiers):
                    qualifiers["zone"] = self._filter_zones(
                        zone_ids=qualifiers["zone"], keep=locations.faces
                    )
                elif zone_ids:
                    zone_ids = self._filter_zones(
                        zone_ids=zone_ids, keep=locations.faces
                    )
                else:
                    if not self._model._server.meet_version("7.1"):
                        raise ValueError(
                            "Querying an ElementalAndFaces result on faces "
                            "currently requires the use of face zone ids in the "
                            "'qualifiers' or the 'zone_ids' arguments."
                        )
                    else:
                        # ElementalAndFaces results have been requested on faces
                        # without defining zones
                        # Do nothing unless face_ids is not defined, in which case we set it to all
                        if face_ids is None:
                            face_ids = self.mesh.face_ids

            elif location == locations.elemental:
                # CFF only returns results on face zones if qualifiers have been set
                if qualifiers and ("zone" in qualifiers):
                    qualifiers["zone"] = self._filter_zones(
                        zone_ids=qualifiers["zone"], keep=locations.elemental
                    )
                elif zone_ids:
                    zone_ids = self._filter_zones(
                        zone_ids=zone_ids, keep=locations.elemental
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

        wf, comp, base_name = self._get_result_workflow(
            base_name=base_name,
            location=location,
            category=category,
            components=components,
            norm=norm,
            selection=selection,
            set_ids=set_ids,
            zone_ids=zone_ids,
            phases=phases,
            species=species,
            qualifiers=qualifiers,
        )

        # Evaluate  the workflow
        fc = wf.get_output("out", dpf.types.fields_container)
        # print(fc)
        if location is None and len(fc) > 0:
            location = fc[0].location
        if location == locations.elemental:
            location = "cells"

        _, _, columns = self._create_components(base_name, category, components)
        return self._create_dataframe(
            fc, location, columns, comp, base_name.split("::")[-1], None
        )

    def _try_get_result_info(self, name: str):
        try:
            return self.result_info[name]
        except ValueError:
            raise ValueError(f"Result {name} is not available.")

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
        location: Union[locations, str, None] = None,
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
            List of faces IDs to get results for. For CFX files, these IDs correspond
            to the COMPOSITE that gathers the related face zones.
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("density")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("density")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def density_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids` are mutually
        exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("density")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("density")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("dynamic_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        """Extract dynamic viscosity results on nodes from the simulation.

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
        result_info = self._try_get_result_info("dynamic_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def dynamic_viscosity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("dynamic_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("dynamic_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("enthalpy")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("enthalpy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def enthalpy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("enthalpy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("enthalpy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("entropy")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("entropy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def entropy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("entropy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("entropy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("epsilon")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("epsilon")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def epsilon_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("epsilon")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("epsilon")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("mach_number")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mach_number")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def mach_number_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mach_number")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        """Extract mach number results on cells from the simulation.

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
        result_info = self._try_get_result_info("mach_number")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def mass_flow_rate(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        .. note::
            This is an integrated result only available on faces.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mass_flow_rate")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            integrated=True,
            native_location=result_info.native_location,
        )

    def mass_flow_rate_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        .. note::
            This is an integrated result only available on faces.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mass_flow_rate")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            integrated=True,
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("mass_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mass_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def mass_fraction_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mass_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mass_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("mean_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        """Extract mean static pressure results on nodes from the simulation.

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
        result_info = self._try_get_result_info("mean_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def mean_static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mean_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mean_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("mean_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mean_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def mean_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mean_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mean_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("mean_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mean_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def mean_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("mean_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("mean_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("omega")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("omega")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def omega_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("omega")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("omega")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("rms_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("rms_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def rms_static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("rms_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("rms_static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("rms_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("rms_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def rms_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("rms_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("rms_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("rms_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("rms_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def rms_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("rms_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("rms_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("specific_heat")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("specific_heat")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def specific_heat_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("specific_heat")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("specific_heat")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def static_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("static_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("superficial_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("superficial_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def superficial_velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("superficial_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("superficial_velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def surface_heat_rate(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        .. note::
            This is an integrated result only available on faces.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, `face_ids`, and `node_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("surface_heat_rate")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            integrated=True,
            native_location=result_info.native_location,
        )

    def surface_heat_rate_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        .. note::
            This is an integrated result only available on faces.

        Arguments `selection`, `set_ids`, `all_sets`, and `times` are mutually
        exclusive.
        If none of the above is given, only the last result will be returned.

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("surface_heat_rate")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            integrated=True,
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("temperature")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("thermal_conductivity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("thermal_conductivity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def thermal_conductivity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("thermal_conductivity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("thermal_conductivity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("total_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("total_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def total_pressure_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("total_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("total_pressure")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("total_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("total_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def total_temperature_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("total_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("total_temperature")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("turbulent_kinetic_energy")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=self.result_info[
                "turbulent_kinetic_energy"
            ].native_location,
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
        result_info = self._try_get_result_info("turbulent_kinetic_energy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=self.result_info[
                "turbulent_kinetic_energy"
            ].native_location,
        )

    def turbulent_kinetic_energy_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("turbulent_kinetic_energy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=self.result_info[
                "turbulent_kinetic_energy"
            ].native_location,
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
        result_info = self._try_get_result_info("turbulent_kinetic_energy")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=self.result_info[
                "turbulent_kinetic_energy"
            ].native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("turbulent_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("turbulent_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def turbulent_viscosity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("turbulent_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("turbulent_viscosity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("velocity")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def velocity_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("velocity")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("volume_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("volume_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def volume_fraction_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("volume_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("volume_fraction")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("wall_shear_stress")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("wall_shear_stress")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def wall_shear_stress_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("wall_shear_stress")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
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
        location: Union[locations, str, None] = None,
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
        location:
            Location to extract results at. Available locations are listed in
            class:`post.locations` and are: `post.locations.nodal`,
            `post.locations.elemental`, and `post.locations.faces`.
            If no location is given, the result is returned as it is stored in the result file.
            Using `post.locations.elemental` gives results with one value for each cell,
            while using `post.locations.nodal` gives results with one value for each node.

        Returns
        -------
        Returns a :class:`ansys.dpf.post.data_object.DataFrame` instance.

        """
        result_info = self._try_get_result_info("y_plus")
        return self._get_result(
            base_name=result_info.operator_name,
            location=location,
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
            native_location=result_info.native_location,
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
        result_info = self._try_get_result_info("y_plus")
        return self._get_result(
            base_name=result_info.operator_name,
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
            native_location=result_info.native_location,
        )

    def y_plus_on_faces(
        self,
        face_ids: Union[List[int], None] = None,
        # cell_ids: Union[List[int], None] = None,
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

        Arguments `selection`, `named_selections`, and `face_ids`
        are mutually exclusive.
        If none of the above is given, results will be extracted for the whole mesh.

        Argument `qualifiers` overrides arguments `zones_ids`, `phases`, and `species`.

        Parameters
        ----------
        face_ids:
            List of IDs of faces to get results for.
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
        result_info = self._try_get_result_info("y_plus")
        return self._get_result(
            base_name=result_info.operator_name,
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
            # cell_ids=cell_ids,
            zone_ids=zone_ids,
            qualifiers=qualifiers,
            phases=phases,
            species=species,
            named_selections=named_selections,
            native_location=result_info.native_location,
        )
