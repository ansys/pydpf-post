"""Module containing the ``Simulation`` class.

Simulation
----------

"""
from abc import ABC
from enum import Enum
from os import PathLike
import re
from typing import Dict, List, Tuple, Union
import warnings

import ansys.dpf.core as dpf
from ansys.dpf.core import DataSources, Model, TimeFreqSupport, Workflow, errors
from ansys.dpf.core.available_result import _result_properties
from ansys.dpf.core.common import elemental_properties
from ansys.dpf.core.plotter import DpfPlotter
from ansys.dpf.core.server_types import BaseServer
import numpy as np

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
from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.meshes import Meshes
from ansys.dpf.post.selection import Selection, _WfNames

component_label_to_index = {
    "1": 0,
    "2": 1,
    "3": 2,
    "4": 3,
    "5": 4,
    "6": 5,
    "X": 0,
    "Y": 1,
    "Z": 2,
    "XX": 0,
    "YY": 1,
    "ZZ": 2,
    "XY": 3,
    "YZ": 4,
    "XZ": 5,
}

vector_component_names = ["X", "Y", "Z"]
matrix_component_names = ["XX", "YY", "ZZ", "XY", "YZ", "XZ"]
principal_names = ["1", "2", "3"]


class ResultCategory(Enum):
    """Enum for available result categories."""

    scalar = 1
    vector = 2
    matrix = 3
    principal = 4
    equivalent = 5


class Simulation(ABC):
    """Base class of all PyDPF-Post simulation types."""

    _vector_component_names = vector_component_names
    _matrix_component_names = matrix_component_names
    _principal_names = principal_names

    def __init__(self, data_sources: DataSources, model: Model):
        """Initialize the simulation using a ``dpf.core.Model`` object."""
        self._model = model
        self._data_sources = data_sources
        self._geometries = []
        self._boundary_conditions = []
        self._loads = []
        self._active_selection = None
        self._named_selections = None
        self._mesh = None
        self._units = None
        self._time_freq_precision = None

    def release_streams(self):
        """Release the streams to data files if any is active."""
        self._model.metadata.release_streams()

    @property
    def results(self) -> List[dpf.result_info.available_result.AvailableResult]:
        r"""Available results.

        Returns a list of available results.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.results) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [...]
        """
        return self._model.metadata.result_info.available_results

    @property
    def result_info(self):
        """Return information concerning the available results."""
        return self._model.metadata.result_info

    @property
    def geometries(self):
        """List of constructed geometries in the simulation.

        Returns a list of geometry objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.geometries) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._geometries

    # @property
    # def boundary_conditions(self):
    #     """List of boundary conditions in the simulation.
    #
    #     Returns a list of boundary_condition objects.
    #
    #     Examples
    #     --------
    #     >>> from ansys.dpf import post
    #     >>> from ansys.dpf.post import examples
    #     >>> simulation = post.load_simulation(examples.static_rst)
    #     >>> print(simulation.boundary_conditions) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    #     []
    #     """
    #     return self._boundary_conditions

    # @property
    # def loads(self):
    #     """List of loads in the simulation.
    #
    #     Returns a list of load objects.
    #
    #     Examples
    #     --------
    #     >>> from ansys.dpf import post
    #     >>> from ansys.dpf.post import examples
    #     >>> simulation = post.load_simulation(examples.static_rst)
    #     >>> print(simulation.loads) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    #     []
    #     """
    #     return self._loads

    @property
    def mesh(self) -> Mesh:
        """Mesh representation of the model.

        Returns a :class:`ansys.dpf.post.mesh.Mesh` object.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
            DPF  Mesh:
              81 nodes
              8 elements
              Unit: m
              With solid (3D) elements
        """
        if self._mesh is None:
            self._mesh = Mesh(self._model.metadata.meshed_region)
        return self._mesh

    @property
    def named_selections(self) -> List[str]:
        """List of named selections in the simulation.

        Returns a list of named selections names.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.named_selections) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['_FIXEDSU']
        """
        if self._named_selections is None:
            self._named_selections = self._model.metadata.available_named_selections
        return self._named_selections

    def plot(
        self,
        mesh: bool = True,
        constructed_geometries: bool = True,
        loads: bool = True,
        boundary_conditions: bool = True,
        **kwargs,
    ):
        """General plot of the simulation object.

        Plots by default the complete mesh contained in the simulation,
        as well as a representation of the constructed geometry,
        the loads, and the boundary conditions currently defined.
        Each representation can be deactivated with its respective boolean argument.

        Parameters
        ----------
        mesh:
            Whether to plot the mesh representation.
        constructed_geometries:
            Whether to plot the constructed geometries.
        loads:
            Whether to plot the loads.
        boundary_conditions:
            Whether to plot the boundary conditions.
        **kwargs
            Additional keyword arguments for the plotter. More information
            are available at :func:`pyvista.plot`.

        Returns
        -------
            Returns a plotter instance of the active visualization backend.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> simulation.plot() # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        """
        plt = DpfPlotter()
        if mesh:
            plt.add_mesh(self.mesh._meshed_region, **kwargs)
        if constructed_geometries:
            for geom in self.geometries:
                getattr(plt, "add_" + str(type(geom).__name__).lower())(geom, **kwargs)
        if loads:
            pass
        if boundary_conditions:
            pass
        plt.show_figure(**kwargs)

    @property
    def active_selection(self) -> Union[Selection, None]:
        """Active selection used by default for result queries.

        Returns a :object:`ansys.dpf.post.selection.Selection` object.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.active_selection = selection
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        """
        return self._active_selection

    @active_selection.setter
    def active_selection(self, selection: Union[Selection, None]):
        """Sets a selection as active on the simulation.

        Activating a given selection on a simulation means it is used
        as a default selection/filter in further result queries.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.active_selection = selection
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        """
        if selection is not None:
            self._active_selection = selection
        else:
            self.deactivate_selection()

    def deactivate_selection(self):
        """Deactivate the currently active selection.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.active_selection = selection
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        >>> simulation.deactivate_selection()
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        None
        """
        self._active_selection = None

    @property
    def _time_frequencies(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._model.metadata.time_freq_support

    @property
    def _get_time_freq_precision(self):
        """Computes a precision for times/frequencies requests based on the underlying support."""
        if self._time_freq_precision is None:
            available_values = self.time_freq_support.time_frequencies.data
            diff = np.diff(available_values, prepend=available_values[0] - 1.0)
            minimum = np.min(diff)
            self._time_freq_precision = minimum / 20.0
        return self._time_freq_precision

    @property
    def time_freq_support(self) -> TimeFreqSupport:
        """Description of the temporal/frequency analysis of the model."""
        return self._time_frequencies

    @property
    def set_ids(self) -> List:
        """Returns the list of set IDs available in the simulation."""
        return list(range(1, self.time_freq_support.n_sets + 1))

    @property
    def units(self):
        """Returns the current units used."""
        us = self._model.metadata.result_info.unit_system
        result_units = us.split(": ")[1].split(", ")
        if self._units is None:
            self._units = {
                "length": result_units[0],
                "mass": result_units[1],
                "force": result_units[2],
                "time": result_units[3],
                "current": result_units[-2],
                "temperature": result_units[-1],
            }
            if len(result_units) > 6:
                self._units["potential"] = result_units[4]
        return self._units

    def split_mesh_by_properties(
        self,
        properties: Union[
            List[elemental_properties],
            Dict[elemental_properties, Union[int, List[int]]],
        ],
    ) -> Meshes:
        """Splits the simulation Mesh according to properties and returns it as Meshes.

        Parameters
        ----------
        properties:
            Elemental properties to split the global mesh by. Returns all meshes if a list,
            or returns only meshes for certain elemental property values if a dict with
            elemental properties labels with associated value or list of values.

        Returns
        -------
        A Meshes entity with resulting meshes.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post.common import elemental_properties
        >>> from ansys.dpf.post import examples
        >>> example_path = examples.download_all_kinds_of_complexity()
        >>> simulation = post.StaticMechanicalSimulation(example_path)
        >>> # Split by elemental properties and get all resulting meshes
        >>> meshes_split = simulation.split_mesh_by_properties(
        ...     properties=[elemental_properties.material,
        ...                 elemental_properties.element_shape]
        ... )
        >>> # Split by elemental properties and only get meshes for certain property values
        >>> # Here: split by material and shape, return only for material 1 and shapes 0 and 1
        >>> meshes_filtered = simulation.split_mesh_by_properties(
        ...     properties={elemental_properties.material: 1,
        ...                 elemental_properties.element_shape: [0, 1]}
        ... )
        """
        split_op = dpf.operators.scoping.split_on_property_type(
            mesh=self._model.metadata.mesh_provider.outputs.mesh,
            requested_location=dpf.locations.elemental,
        )
        values = None
        if isinstance(properties, dict):
            values = properties.values()
            properties = list(properties.keys())
        for i, prop in enumerate(properties):
            split_op.connect(13 + i, prop)
        scopings_container = split_op.outputs.mesh_scoping()
        # meshes = []
        meshes_container = dpf.MeshesContainer()
        for label in scopings_container.labels:
            meshes_container.add_label(label)
        for i, scoping in enumerate(scopings_container):
            mesh_from_scoping = dpf.operators.mesh.from_scoping(
                scoping=scoping,
                mesh=self._model.metadata.mesh_provider.outputs.mesh,
            )
            # meshes.append(mesh_from_scoping.outputs.mesh())
            meshes_container.add_mesh(
                mesh=mesh_from_scoping.outputs.mesh(),
                label_space=scopings_container.get_label_space(i),
            )
        meshes = Meshes(meshes_container=meshes_container)
        if values is None:
            return meshes
        else:
            return meshes.select(**dict(zip(properties, values)))

    def __str__(self):
        """Get the string representation of this class."""
        txt = (
            "%s." % re.sub(r"(?<!^)(?=[A-Z])", " ", type(self).__name__)
            + "\n\n\nData Sources\n------------------------------\n"
        )
        ds_str = self._model._data_sources.__str__()
        txt += ds_str
        txt += "\n\n"
        txt += self._model.__str__()
        return txt

    def _build_components_for_vector(self, base_name, components):
        out, columns = self._build_components(
            base_name, components, self._vector_component_names
        )
        return out, columns

    def _build_components_for_matrix(self, base_name, components):
        out, columns = self._build_components(
            base_name, components, self._matrix_component_names
        )
        return out, columns

    def _build_components(self, base_name, components, component_names):
        # Create operator internal names based on components
        out = []
        if components is None:
            out = None
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if not isinstance(components, list):
                raise ValueError(
                    "Argument 'components' must be an int, a str, or a list of either."
                )
            for comp in components:
                if not (isinstance(comp, str) or isinstance(comp, int)):
                    raise ValueError(
                        "Argument 'components' can only contain integers and/or strings.\n"
                        f"The provided component '{comp}' is not valid."
                    )
                if isinstance(comp, int):
                    comp = str(comp)
                if comp not in component_label_to_index.keys():
                    raise ValueError(
                        f"Component {comp} is not valid. Please use one of: "
                        f"{list(component_label_to_index.keys())}."
                    )
                out.append(component_label_to_index[comp])

        # Take unique values and build names list
        if out is None:
            columns = [base_name + comp for comp in component_names]
        else:
            out = list(set(out))
            columns = [base_name + component_names[i] for i in out]
        return out, columns

    def _build_components_for_principal(self, base_name, components):
        # Create operator internal names based on principal components
        out = []
        if components is None:
            components = [1]

        if isinstance(components, int) or isinstance(components, str):
            components = [components]
        if not isinstance(components, list):
            raise ValueError(
                "Argument 'components' must be an int, a str, or a list of either."
            )
        for comp in components:
            if not (isinstance(comp, str) or isinstance(comp, int)):
                raise ValueError(
                    "Argument 'components' can only contain integers and/or strings."
                )
            if str(comp) not in self._principal_names:
                raise ValueError(
                    "A principal component ID must be one of: "
                    f"{self._principal_names}."
                )
            out.append(int(comp) - 1)

        # Take unique values
        if out is not None:
            out = list(set(out))
        # Build columns names
        if out is None:
            columns = [base_name + str(comp) for comp in self._principal_names]
        else:
            columns = [base_name + self._principal_names[i] for i in out]
        return out, columns

    def _build_result_operator(
        self,
        name: str,
        location: Union[locations, str],
        force_elemental_nodal: bool,
    ) -> dpf.Operator:
        op = self._model.operator(name=name)
        op.connect(7, self.mesh._meshed_region)
        if force_elemental_nodal:
            op.connect(9, "ElementalNodal")
        else:
            op.connect(9, location)
        return op

    def _build_result_workflow(
        self,
        name: str,
        location: Union[locations, str],
        force_elemental_nodal: bool,
    ) -> (dpf.Workflow, dpf.Operator):
        op = self._model.operator(name=name)
        op.connect(7, self.mesh._meshed_region)
        if force_elemental_nodal:
            op.connect(9, "ElementalNodal")
        elif location:
            op.connect(9, location)
        wf = Workflow(server=self._model._server)
        wf.add_operator(op)
        wf.set_input_name(_WfNames.read_cyclic, op, 14)
        wf.set_input_name(_WfNames.cyclic_sectors_to_expand, op, 18)
        wf.set_input_name(_WfNames.cyclic_phase, op, 19)
        wf.set_output_name(_WfNames.result, op, 0)
        return wf, op

    def _append_norm(self, wf, out, base_name):
        """Append a norm operator to the current result workflow."""
        norm_op = self._model.operator(name="norm_fc")
        norm_op.connect(0, out)
        wf.add_operator(operator=norm_op)
        base_name += "_N"
        out = norm_op.outputs.fields_container
        comp = None
        return wf, out, comp, base_name

    def _create_components(self, base_name, category, components):
        comp = None
        # Build the list of requested results
        if category in [ResultCategory.scalar, ResultCategory.equivalent]:
            # A scalar or equivalent result has no components
            to_extract = None
            columns = [base_name]
        elif category == ResultCategory.vector:
            # A vector result can have components selected
            to_extract, columns = self._build_components_for_vector(
                base_name=base_name, components=components
            )
            if to_extract is not None:
                comp = [self._vector_component_names[i] for i in to_extract]
            else:
                comp = self._vector_component_names
        elif category == ResultCategory.matrix:
            # A vector result can have components selected
            to_extract, columns = self._build_components_for_matrix(
                base_name=base_name, components=components
            )
            if to_extract is not None:
                comp = [self._matrix_component_names[i] for i in to_extract]
            else:
                comp = self._matrix_component_names
        elif category == ResultCategory.principal:
            # A principal type of result can have components selected
            to_extract, columns = self._build_components_for_principal(
                base_name=base_name, components=components
            )
            comp = [self._principal_names[i] for i in to_extract]
        else:
            raise ValueError(f"'{category}' is not a valid category value.")
        return comp, to_extract, columns

    def _generate_disp_workflow(self, fc, selection) -> Union[dpf.Workflow, None]:
        # Check displacement is an available result
        if not any(
            [
                result.name == "displacement"
                for result in self._model.metadata.result_info.available_results
            ]
        ):
            return None
        # Build an equivalent workflow for displacement for plots and animations
        disp_wf = dpf.Workflow(server=fc._server)

        disp_op = dpf.operators.result.displacement(
            data_sources=self._model.metadata.data_sources,
            streams_container=self._model.metadata.streams_provider,
            server=fc._server,
        )
        # Connect time_scoping (do not connect mesh_scoping as we want to deform the whole mesh)
        disp_wf.set_input_name("time_scoping", disp_op.inputs.time_scoping)
        disp_wf.connect_with(
            selection.time_freq_selection._selection,
            output_input_names=("scoping", "time_scoping"),
        )

        # Shell layer selection step
        shell_layer_op = dpf.operators.utility.change_shell_layers(
            fields_container=disp_op.outputs.fields_container,
            server=fc._server,
        )
        # Expose shell layer input as workflow input
        disp_wf.set_input_name("shell_layer_int", shell_layer_op.inputs.e_shell_layer)

        # Merge shell and solid fields at same set
        merge_op = dpf.operators.utility.merge_fields_by_label(
            fields_container=shell_layer_op.outputs.fields_container_as_fields_container,
            label="eltype",
            server=fc._server,
        )

        # Expose output
        disp_wf.set_output_name("output", merge_op.outputs.fields_container)

        return disp_wf

    def _create_dataframe(
        self, fc, location, columns, comp, base_name, disp_wf=None, submesh=None
    ):
        # Test for empty results
        if (len(fc) == 0) or all([len(f) == 0 for f in fc]):
            warnings.warn(
                message=f"Returned Dataframe with columns {columns} is empty.",
                category=UserWarning,
            )
        unit = None
        times = [""]
        if len(fc) > 0:
            unit = fc[0].unit
            times = fc.get_available_ids_for_label("time")
            if submesh is not None:
                for i_field in range(len(fc)):
                    bind_support_op = dpf.operators.utility.bind_support(
                        fc[i_field],
                        submesh,
                        server=fc._server,
                    )
                    fc.add_field(fc.get_label_space(i_field), bind_support_op.eval())
        if unit == "":
            unit = None
        comp_index = None
        if comp is not None:
            comp_index = CompIndex(values=comp)
        row_indexes = [MeshIndex(location=location, fc=fc)]
        if comp_index is not None:
            row_indexes.append(comp_index)
        if len(fc) > 0:
            times = fc.get_available_ids_for_label("time")
        column_indexes = [ResultsIndex(values=[base_name], units=[unit])]
        if times:
            column_indexes.append(SetIndex(values=times))
        label_indexes = []
        # Get the label name values for each label
        for label in fc.labels:
            # Do not treat time
            if label in ["time"]:
                continue
            if len(fc) > 0:
                # Get ID values for this label
                values = fc.get_available_ids_for_label(label)
                # Then try to gather the correspond string values for display
                try:
                    label_support = self.result_info.qualifier_label_support(label)
                    names_field = label_support.string_field_support_by_property(
                        "names"
                    )
                    values = [
                        names_field.get_entity_data_by_id(value)[0] + f" ({value})"
                        for value in values
                    ]
                except (
                    ValueError,
                    errors.DPFServerException,
                    errors.DpfVersionNotSupported,
                    AttributeError,  # Error in GitHub doc CI (Linux) where label_support = None
                ):
                    values = None
            else:
                values = []
            label_indexes.append(LabelIndex(name=label, values=values))

        column_indexes.extend(label_indexes)
        column_index = MultiIndex(indexes=column_indexes)

        row_index = MultiIndex(
            indexes=row_indexes,
        )

        df = DataFrame(
            data=fc,
            columns=column_index,
            index=row_index,
        )
        df._disp_wf = disp_wf

        # Return the result wrapped in a DPF_Dataframe
        return df

    @staticmethod
    def _treat_cyclic(expand_cyclic, phase_angle_cyclic, result_wf):
        if expand_cyclic is not False:
            # If expand_cyclic is a list
            if isinstance(expand_cyclic, list) and len(expand_cyclic) > 0:
                # If a list of sector numbers, directly connect it to the num_sectors pin
                if all(
                    [
                        isinstance(expand_cyclic_i, int)
                        for expand_cyclic_i in expand_cyclic
                    ]
                ):
                    if any([i < 1 for i in expand_cyclic]):
                        raise ValueError(
                            "Sector selection with 'expand_cyclic' starts at 1."
                        )
                    result_wf.connect(
                        _WfNames.cyclic_sectors_to_expand,
                        [i - 1 for i in expand_cyclic],
                    )
                # If any is a list, treat it as per stage num_sectors
                elif any(
                    [
                        isinstance(expand_cyclic_i, list)
                        for expand_cyclic_i in expand_cyclic
                    ]
                ):
                    # Create a ScopingsContainer to fill
                    sectors_scopings = dpf.ScopingsContainer()
                    sectors_scopings.labels = ["stage"]
                    # For each potential num_sectors, check either an int or a list of ints
                    for i, num_sectors_stage_i in enumerate(expand_cyclic):
                        # Prepare num_sectors data
                        if isinstance(num_sectors_stage_i, int):
                            num_sectors_stage_i = [num_sectors_stage_i]
                        elif isinstance(num_sectors_stage_i, list):
                            if not all(
                                [isinstance(n, int) for n in num_sectors_stage_i]
                            ):
                                raise ValueError(
                                    "'expand_cyclic' only accepts lists of int values >= 1."
                                )
                        # num_sectors_stage_i is now a list of int,
                        # add an equivalent Scoping with the correct 'stage' label value
                        if any([i < 1 for i in num_sectors_stage_i]):
                            raise ValueError(
                                "Sector selection with 'expand_cyclic' starts at 1."
                            )
                        sectors_scopings.add_scoping(
                            {"stage": i},
                            dpf.Scoping(ids=[i - 1 for i in num_sectors_stage_i]),
                        )
                    result_wf.connect(
                        _WfNames.cyclic_sectors_to_expand, inpt=sectors_scopings
                    )
            elif not isinstance(expand_cyclic, bool):
                raise ValueError(
                    "'expand_cyclic' argument can only be a boolean or a list."
                )
            result_wf.connect(_WfNames.read_cyclic, 3)  # Connect the read_cyclic pin
        else:
            result_wf.connect(_WfNames.read_cyclic, 1)  # Connect the read_cyclic pin
        if phase_angle_cyclic is not None:
            if isinstance(phase_angle_cyclic, int):
                phase_angle_cyclic = float(phase_angle_cyclic)
            if not isinstance(phase_angle_cyclic, float):
                raise ValueError(
                    "'phase_angle_cyclic' argument only accepts a single float value."
                )
            result_wf.connect(_WfNames.cyclic_phase, phase_angle_cyclic)
        return result_wf


class MechanicalSimulation(Simulation, ABC):
    """Base class for mechanical type simulations.

    This class provides common methods and properties for all mechanical type simulations.
    """

    def __init__(
        self,
        result_file: Union[PathLike, str, DataSources],
        server: Union[BaseServer, None] = None,
    ):
        """Instantiate a mechanical type simulation."""
        model = dpf.Model(result_file, server=server)
        data_sources = model.metadata.data_sources
        super().__init__(data_sources=data_sources, model=model)

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
        external_layer: bool = False,
        skin: Union[bool, List[int]] = False,
        expand_cyclic: Union[bool, List[Union[int, List[int]]]] = True,
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
        tot = (skin is not None and skin is not False) + (
            external_layer is not None and external_layer is not False
        )
        if tot > 1:
            raise ValueError(
                "Arguments selection, skin, and external_layer are mutually exclusive"
            )
        if selection is not None:
            return selection
        else:
            selection = Selection(server=self._model._server)
        # Create the SpatialSelection

        # First: the skin and the external layer to be able to have both a mesh scoping and
        # the skin/external layer
        if (skin is not None and skin is not False) or (
            external_layer is not None and external_layer is not False
        ):
            res = (
                _result_properties[base_name]
                if base_name in _result_properties
                else None
            )
            location = (
                locations.elemental_nodal
                if self._requires_manual_averaging(base_name, location, category, None)
                else location
            )
            if external_layer not in [None, False]:
                selection.select_external_layer(
                    elements=external_layer if external_layer is not True else None,
                    location=location,
                    result_native_location=res["location"]
                    if res is not None
                    else location,
                    is_model_cyclic=self._model.operator("is_cyclic").eval()
                    if expand_cyclic is not False
                    else "not_cyclic",
                )
            else:
                selection.select_skin(
                    elements=skin if skin is not True else None,
                    location=location,
                    result_native_location=res["location"]
                    if res is not None
                    else location,
                    is_model_cyclic=self._model.operator("is_cyclic").eval()
                    if expand_cyclic is not False
                    else "not_cyclic",
                )
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
                    raise ValueError(f"Could not find time={t} in the simulation.")
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

    def _requires_manual_averaging(
        self,
        base_name: str,
        location: str,
        category: ResultCategory,
        selection: Selection,
    ):
        res = _result_properties[base_name] if base_name in _result_properties else None
        if category == ResultCategory.equivalent and base_name[0] == "E":  # strain eqv
            return True
        if res is not None and selection is not None:
            return selection.requires_manual_averaging(
                location=location,
                result_native_location=res["location"],
                is_model_cyclic=self._model.operator("is_cyclic").eval(),
            )
        return False

    def _create_averaging_operator(
        self,
        location: str,
        selection: Selection,
    ):
        average_op = None
        first_average_op = None
        forward = None
        if _WfNames.skin in selection.spatial_selection._selection.output_names:
            if self._model._server.meet_version("6.2"):
                first_average_op = self._model.operator(name="solid_to_skin_fc")
                forward = first_average_op
            else:
                first_average_op = self._model.operator(name="solid_to_skin")
                forward = self._model.operator(name="forward_fc")
                forward.connect(0, first_average_op, 0)
            average_wf = dpf.Workflow(server=self._model._server)
            if hasattr(first_average_op.inputs, "mesh_scoping"):
                inpt = (
                    first_average_op.inputs.mesh_scoping
                )  # To keep for retro-compatibility
            else:
                inpt = first_average_op.inputs.mesh
            average_wf.set_input_name(_WfNames.skin, inpt)
            average_wf.connect_with(
                selection.spatial_selection._selection,
                output_input_names={_WfNames.skin: _WfNames.skin},
            )

        if location == locations.nodal:
            average_op = self._model.operator(name="to_nodal_fc")
        elif location == locations.elemental:
            average_op = self._model.operator(name="to_elemental_fc")
        if average_op and forward:
            average_op.connect(0, forward, 0)
        else:
            first_average_op = average_op

        if first_average_op is not None and average_op is not None:
            return (first_average_op, average_op)
