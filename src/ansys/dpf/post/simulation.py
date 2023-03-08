"""Module containing the ``Simulation`` class.

Simulation
----------

"""
from abc import ABC
from enum import Enum
from os import PathLike
import re
from typing import List, Tuple, Union
import warnings

import ansys.dpf.core as dpf
from ansys.dpf.core import DataSources, Model, TimeFreqSupport
from ansys.dpf.core.plotter import DpfPlotter
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
from ansys.dpf.post.selection import Selection

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
        self._units = {
            "time/frequency": self.time_freq_support.time_frequencies.unit,
            "distance": self._model.metadata.meshed_region.unit,
        }
        self._time_freq_precision = None

    def release_streams(self):
        """Release the streams to data files if any is active."""
        self._model.metadata.release_streams()

    @property
    def results(self) -> List[str]:
        r"""Available results.

        Returns a list of available results as strings.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.results) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['displacement\nOperator name: "U"\n...Units: degc\n']
        """
        return [
            str(result) for result in self._model.metadata.result_info.available_results
        ]

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
            plt.add_mesh(self.mesh._meshed_region)
        if constructed_geometries:
            for geom in self.geometries:
                getattr(plt, "add_" + str(type(geom).__name__).lower())(geom)
        if loads:
            pass
        if boundary_conditions:
            pass
        plt.show_figure()

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
        """Returns the current time/frequency and distance units used."""
        return self._units

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
        )

        # Expose output
        disp_wf.set_output_name("output", merge_op.outputs.fields_container)

        return disp_wf

    def _create_dataframe(self, fc, location, columns, comp, base_name, disp_wf=None):
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

        df = DataFrame(
            data=fc,
            columns=column_index,
            index=row_index,
        )
        df._disp_wf = disp_wf

        # Return the result wrapped in a DPF_Dataframe
        return df


class MechanicalSimulation(Simulation, ABC):
    """Base class for mechanical type simulations.

    This class provides common methods and properties for all mechanical type simulations.
    """

    def __init__(self, result_file: Union[PathLike, str]):
        """Instantiate a mechanical type simulation."""
        model = dpf.Model(result_file)
        data_sources = model.metadata.data_sources
        super().__init__(data_sources=data_sources, model=model)

    def _build_selection(
        self,
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
