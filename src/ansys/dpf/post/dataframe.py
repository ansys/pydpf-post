"""Module containing the ``DataFrame`` class."""
from __future__ import annotations

import itertools
from os import PathLike
from typing import List, Union
import warnings

import ansys.dpf.core as dpf
from ansys.dpf.core.common import shell_layers
from ansys.dpf.core.dpf_array import DPFArray

from ansys.dpf.post import locations
from ansys.dpf.post.index import (
    CompIndex,
    Index,
    LabelIndex,
    MeshIndex,
    MultiIndex,
    ResultsIndex,
    SetIndex,
)

display_width = 80
display_max_colwidth = 10
display_max_lines = 6


class DataFrame:
    """A DataFrame style API to manipulate DPF data."""

    def __init__(
        self,
        data: dpf.FieldsContainer,
        index: Union[MultiIndex, Index, List[int]],
        columns: Union[MultiIndex, Index, List[str], None] = None,
    ):
        """Creates a DPF DataFrame based on the given input data.

        Parameters
        ----------
        data:
            Data to use.
        index:
            Row indexing (labels) to use.
        columns:
            Column indexing (labels) to use.
        """
        self._index = index
        if isinstance(data, dpf.FieldsContainer):
            self._fc = data
            # if index is None:
            #     raise NotImplementedError("Creation from FieldsContainer without index "
            #                               "is not yet supported")
            #     # self._index = Index(
            #     #     name=location_to_label[data[0].location], values=None
            #     # )
        else:
            raise ValueError(
                f"Input data type '{type(data)}' is not a valid data type for DataFrame creation."
            )
        if columns is not None:
            self._columns = columns
        else:
            self._columns = None

        # if parent_simulation is not None:
        #     self._parent_simulation = weakref.ref(parent_simulation)

        self._str = None
        self._last_display_width = display_width
        self._last_display_max_colwidth = display_max_colwidth

    @property
    def columns(self):
        """Returns the column labels of the DataFrame."""
        if self._columns is None:
            indexes = [ResultsIndex(values=[self._fc[0].name.split("_")])]
            indexes.extend(
                [
                    LabelIndex(name=label, values=self._fc.get_label_scoping(label).ids)
                    for label in self._fc.labels
                ]
            )
            self._columns = MultiIndex(indexes=indexes)
        return self._columns

    @property
    def index(self) -> Union[MultiIndex, Index]:
        """Returns the Index or MultiIndex for the rows of the DataFrame."""
        return self._index

    @property
    def axes(self) -> List[MultiIndex]:
        """Returns a list with the row MultiIndex first and the columns MultiIndex second."""
        return [self.index, self.columns]

    @property
    def results_index(self) -> Union[ResultsIndex, None]:
        """Returns the available ResultsIndex is present."""
        results_index = self.columns.results_index
        if results_index is None:
            results_index = self.index.results_index
        return results_index

    @property
    def mesh_index(self) -> Union[MeshIndex, None]:
        """Returns the available MeshIndex is present."""
        mesh_index = self.columns.mesh_index
        if mesh_index is None:
            mesh_index = self.index.mesh_index
        return mesh_index

    @property
    def labels(self) -> List[str]:
        """Returns a list of the names of available LabelIndex indexes."""
        names = [index.name for index in self.index if isinstance(index, LabelIndex)]
        names.extend(
            [index.name for index in self.columns if isinstance(index, LabelIndex)]
        )
        return names

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.FieldsContainer` object."""
        return self._fc

    def _validate_arguments(self, arguments):
        """Check for invalid arguments based on available Index names."""
        rows, columns = self.axes
        axes_names = rows.names
        axes_names.extend(columns.names)
        for argument in arguments.keys():
            if argument not in axes_names:
                raise ValueError(
                    f"The DataFrame has no axis {argument}, cannot select it. "
                    f"Available axes are: {axes_names}."
                )

    def select(self, **kwargs) -> DataFrame:
        """Returns a new DataFrame based on selection criteria (value-based).

        Parameters
        ----------
        **kwargs:
            This function accepts as argument any of the Index names available associated with a
            value or a list of values.
            For example, if 'time' is an available class:`Index <ansys.dpf.post.index.Index>`
            of the class:`DataFrame <ansys.dpf.post.DataFrame>` `df`, then you can select the time 1
             by using `df.select(time=1)`.
            One can get the list of available axes using
            :func:`DataFrame.axes <ansys.dpf.post.DataFrame.axes>`.

        Returns
        -------
            A DataFrame of the selected values.

        """
        self._validate_arguments(arguments=kwargs)
        if "set_id" in kwargs.keys():
            kwargs["time"] = kwargs["set_id"]
        # Initiate a workflow
        server = self._fc._server
        wf = dpf.Workflow(server=server)
        wf.progress_bar = False
        input_fc = self._fc
        fc = self._fc
        out = None

        # Initiate future DataFrame indexes with those from initial DataFrame
        mesh_index = None
        comp_index = None
        for index in self.index:
            if isinstance(index, MeshIndex):
                mesh_index = index
            if isinstance(index, CompIndex):
                comp_index = index
        location = self._fc[0].location
        results_index = self.results_index

        # Treat selection on a label
        if any([label in kwargs.keys() for label in self._fc.labels]):
            fc = dpf.FieldsContainer()
            fc.labels = self._fc.labels
            for i, field in enumerate(self._fc):
                label_space = self._fc.get_label_space(i)
                for fc_key in label_space.keys():
                    if fc_key in kwargs.keys() or (
                        fc_key == "time" and "set_id" in kwargs.keys()
                    ):
                        key = fc_key
                        if fc_key == "time" and "set_id" in kwargs.keys():
                            key = "set_id"
                        if not isinstance(kwargs[key], list):
                            kwargs[key] = [kwargs[key]]
                        if label_space[fc_key] not in kwargs[key]:
                            continue
                    fc.add_field(label_space=label_space, field=field)
            input_fc = fc

        # # Treat selection on components
        if "comp" in kwargs.keys():
            from ansys.dpf.post.simulation import component_label_to_index

            comp_to_extract = kwargs["comp"]
            if not isinstance(comp_to_extract, list):
                comp_to_extract = [comp_to_extract]
            component_indexes = [component_label_to_index[c] for c in comp_to_extract]
            selector_fc = dpf.operators.logic.component_selector_fc(
                fields_container=input_fc,
                component_number=component_indexes,
                server=server,
            )
            out = selector_fc.outputs.fields_container
            comp_index = CompIndex(values=comp_to_extract)

        # Treat selection on DataFrame.index (rescope)
        if isinstance(self.index, MeshIndex):
            mesh_index_name = self.index.name
        else:
            mesh_index_name = self.index.mesh_index.name
        if (
            mesh_index_name in kwargs.keys()
            and mesh_index.location != locations.elemental_nodal
        ):
            if "node" in mesh_index_name:
                node_ids = kwargs[mesh_index_name]
                if not isinstance(node_ids, (DPFArray, list)):
                    node_ids = [node_ids]
                mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
                    node_ids=node_ids,
                    server=server,
                )
            elif "element" in mesh_index_name:
                element_ids = kwargs[mesh_index_name]
                if not isinstance(element_ids, list):
                    element_ids = [element_ids]
                mesh_scoping = dpf.mesh_scoping_factory.elemental_scoping(
                    element_ids=element_ids,
                    server=server,
                )
            else:
                raise NotImplementedError(
                    f"Selection on a DataFrame with index "
                    f"'{mesh_index_name}' is not yet supported"
                )
            rescope_fc = dpf.operators.scoping.rescope_fc(
                fields_container=input_fc,
                mesh_scoping=mesh_scoping,
                server=server,
            )
            out = rescope_fc.outputs.fields_container
            mesh_index = MeshIndex(location=location, values=mesh_scoping.ids)
        elif (
            mesh_index_name in kwargs.keys()
            and mesh_index.location == locations.elemental_nodal
        ):
            raise NotImplementedError(
                f"Element selection on a DataFrame with elemental nodal results "
                f"is not yet supported"
            )

        if out is not None:
            wf.set_output_name("out", out)
            fc = wf.get_output("out", dpf.FieldsContainer)

        row_indexes = [mesh_index]
        if comp_index is not None:
            row_indexes.append(comp_index)
        row_index = MultiIndex(
            indexes=row_indexes,
        )

        column_indexes = [
            results_index,
            SetIndex(values=fc.get_available_ids_for_label("time")),
        ]
        label_indexes = []
        for label in fc.labels:
            if label not in ["time"]:
                column_indexes.append(
                    LabelIndex(name=label, values=fc.get_available_ids_for_label(label))
                )
        column_indexes.extend(label_indexes)
        column_index = MultiIndex(indexes=column_indexes)

        return DataFrame(
            data=fc,
            columns=column_index,
            index=row_index,
        )

    def iselect(self, **kwargs) -> DataFrame:
        """Returns a new DataFrame based on selection criteria (index-based).

        Parameters
        ----------
        **kwargs:
            This function accepts as argument any of the Index names available associated with a
            value or a list of values.
            For example, if 'time' is an available class:`Index <ansys.dpf.post.index.Index>`
            of the class:`DataFrame <ansys.dpf.post.DataFrame>` `df`, then you can select the first
            `time` value by using `df.select(set_id=0)`.
            One can get the list of available axes using
            :func:`DataFrame.axes <ansys.dpf.post.DataFrame.axes>`.

        Returns
        -------
            A DataFrame of the selected values.

        """
        self._validate_arguments(arguments=kwargs)
        for label in kwargs.keys():
            indices = kwargs[label]
            if label in self.index.names:
                ids = getattr(self.index, label).values[indices]
            else:
                ids = getattr(self.columns, label).values[indices]
            if isinstance(ids, DPFArray):
                ids = ids.tolist()
            kwargs[label] = ids
        return self.select(**kwargs)

    def __len__(self):
        """Return the length of the DataFrame."""
        return len(self.index)

    def __str__(self) -> str:
        """String representation of the DataFrame."""
        if (
            (self._str is None)
            or (self._last_display_width != display_width)
            or (self._last_display_max_colwidth != display_max_colwidth)
        ):
            self._update_str(width=display_width, max_colwidth=display_max_colwidth)
            self._last_display_width = display_width
            self._last_display_max_colwidth = display_max_colwidth
        return self._str

    def _update_str(self, width: int, max_colwidth: int):
        """Updates the DataFrame string representation using given display options.

        The string representation is limited to five lines, meaning any DataFrame with more than
        five rows will be truncated and only the five first rows are displayed.

        Parameters
        ----------
        width:
            Number of characters to use for the total width.
        max_colwidth:
            Maximum number of characters to use for each column.
        """
        lines = []
        empty = " " * max_colwidth
        truncated_str = "...".rjust(max_colwidth)
        max_n_col = width // max_colwidth
        max_n_rows = display_max_lines
        # Create lines with row labels and values
        num_column_indexes = len(self.columns)
        num_rows_indexes = len(self.index)
        n_max_value_col = max_n_col - num_rows_indexes
        column_headers = []
        for column_index in self.columns:
            column_headers.append(
                empty * (num_rows_indexes - 1) + column_index.name.rjust(max_colwidth)
            )
        lines.extend(column_headers)

        row_headers = "".join(
            [row_index.name.rjust(max_colwidth) for row_index in self.index]
        )
        lines.append(row_headers)
        entity_ids = []
        truncated = False
        num_mesh_entities_to_ask = self._fc[0].size
        if num_mesh_entities_to_ask > max_n_rows:
            num_mesh_entities_to_ask = max_n_rows
            truncated = True
        lists = []
        # Create combinations for rows
        entity_ids = None
        for index in self.index:
            if isinstance(index, MeshIndex):
                if index._values is not None:
                    values = index.values
                else:
                    values = self._first_n_ids_first_field(num_mesh_entities_to_ask)
                entity_ids = values
            elif isinstance(index, CompIndex):
                values = index.values
            else:
                values = index.values
            lists.append(values)
        row_combinations = [p for p in itertools.product(*lists)][:max_n_rows]

        # Add row headers for the first combinations (max_n_rows)
        previous_combination = [None] * len(lists)
        for combination in row_combinations:
            line = "".join(
                [
                    str(combination[i]).rjust(max_colwidth)
                    if combination[i] != previous_combination[i]
                    else empty
                    for i in range(len(combination))
                ]
            )
            previous_combination = combination
            lines.append(line)

        # Create combinations for columns
        entity_ids = None
        lists = []
        label_positions_in_combinations = {}
        comp_values = [1]
        for position, index in enumerate(self.columns):
            if isinstance(index, MeshIndex):
                if index._values is not None:
                    values = index.values
                    entity_ids = values
                else:
                    values = self._first_n_ids_first_field(num_mesh_entities_to_ask)
                    entity_ids = values
            elif isinstance(index, CompIndex):
                values = index.values
                comp_values = values
            else:
                values = index.values
            label_positions_in_combinations[index.name] = position
            lists.append(values)
        combinations = [p for p in itertools.product(*lists)]
        truncated_columns = False
        if len(combinations) > n_max_value_col:
            truncated_columns = True
            combinations = combinations[:n_max_value_col]

        def flatten(arr):
            new_arr = []
            for item in arr:
                if isinstance(item, list):
                    new_arr.extend(flatten(item))
                else:
                    new_arr.append(item)
            return new_arr

        def treat_elemental_nodal(treat_lines, pos, n_comp, n_ent, n_lines):
            # Update row headers
            elem_headers = treat_lines[pos : pos + n_comp]
            new_elem_headers = []
            for i_ent in range(1, n_ent + 1):
                for header in elem_headers:
                    new_elem_header = [
                        header[:max_colwidth]
                        + header[max_colwidth + 4 :]
                        + f" ({i_ent})"
                    ]
                    # elem_header_i.extend(elem_headers[1:])
                    new_elem_headers.extend(new_elem_header)
                treat_lines = treat_lines[:pos] + new_elem_headers + treat_lines[pos:]
            return treat_lines[:n_lines]

        # Add text for the first n_max_value_col columns
        previous_combination = [None] * len(lists)
        for i_c, combination in enumerate(combinations[:n_max_value_col]):
            to_append = [
                str(combination[i]).rjust(max_colwidth)
                if combination[i] != previous_combination[i]
                else empty
                for i in range(len(combination))
            ]
            to_append.append(empty)
            # Get data in the FieldsContainer for those positions
            # Create label_space from combination
            label_space = {}
            for label_name in self.labels:
                value = combination[label_positions_in_combinations[label_name]]
                if label_name == "set_id":
                    label_name = "time"
                label_space[label_name] = value
            fields = self._fc.get_fields(label_space=label_space)
            values = []
            if entity_ids is None:
                for field in fields:
                    array_values = []
                    position = len(column_headers) + 1
                    for k in list(range(num_mesh_entities_to_ask)):
                        try:
                            values_list = field.get_entity_data(k).tolist()
                        except Exception as e:
                            values_list = [[None] * len(comp_values)]
                        num_entities = len(values_list)
                        if isinstance(values_list[0], list):
                            num_components = len(values_list[0])
                        else:
                            num_components = 1
                        current_number_lines = len(lines)
                        # Detect number of nodes when elemental nodal and update headers to
                        # repeat for each node (do not print their ID for now)
                        if (
                            i_c == 0
                            and field.location == locations.elemental_nodal
                            and len(values_list) != 0
                        ):
                            lines = treat_elemental_nodal(
                                lines,
                                position,
                                num_components,
                                num_entities,
                                current_number_lines,
                            )
                        array_values.append(values_list)
                        if (
                            position + num_entities * num_components
                            >= current_number_lines + len(column_headers) + 1
                        ):
                            if k < num_mesh_entities_to_ask:
                                truncated = True
                        if position >= current_number_lines:
                            break
                        position += num_entities * num_components
                    if array_values:
                        array_values = flatten(array_values)
                        values.extend(array_values)
            else:
                array_values = []
                position = len(column_headers) + 1
                for entity_id in entity_ids:
                    for field in fields:
                        try:
                            values_list = field.get_entity_data_by_id(
                                entity_id
                            ).tolist()
                        except Exception as e:
                            values_list = [[None] * len(comp_values)]
                        num_entities = len(values_list)
                        num_components = len(values_list[0])
                        current_number_lines = len(lines)
                        # Detect number of nodes when elemental nodal and update headers to
                        # repeat for each node (do not print their ID for now)
                        if (
                            i_c == 0
                            and field.location == locations.elemental_nodal
                            and len(values_list) != 0
                        ):
                            lines = treat_elemental_nodal(
                                lines,
                                position,
                                num_components,
                                num_entities,
                                current_number_lines,
                            )
                        array_values.append(values_list)
                        if (
                            position + num_entities * num_components
                            >= current_number_lines + len(column_headers) + 1
                        ):
                            if k < num_mesh_entities_to_ask:
                                truncated = True
                        if position >= current_number_lines:
                            break
                        position += num_entities * num_components
                        if array_values:
                            array_values = flatten(array_values)
                            values.extend(array_values)

            # take_comp_map = [True] * len(values)
            # if comp_values is not None:

            value_strings = []
            for value in values[: len(lines) - (num_column_indexes + 1)]:
                if value is not None:
                    value_string = f"{value:.{max_colwidth - 8}e}".rjust(max_colwidth)
                else:
                    value_string = empty
                value_strings.append(value_string)
            to_append.extend(value_strings)
            previous_combination = combination
            # print(to_append)
            # print(len(to_append))
            # print(len(lines))
            for i in range(len(lines)):
                lines[i] = lines[i] + to_append[i]

        if truncated_columns:
            for i in range(len(lines)):
                lines[i] = lines[i] + truncated_str

        if truncated:
            lines.append(truncated_str)
        txt = "\n" + "".join([line + "\n" for line in lines])
        self._str = txt

    def _first_n_ids_first_field(self, n: int):
        """Return the n first entity IDs of the first field in the underlying FieldsContainer."""
        return self._fc[0].scoping.ids[:n]

    def __repr__(self):
        """Representation of the DataFrame."""
        return f"DataFrame<index={self.index}, columns={self.columns}>"

    def plot(self, shell_layer=shell_layers.top, **kwargs):
        """Plot the result.

        Parameters
        ----------
        shell_layer:
            Shell layer to show if multi-layered shell data present. Defaults to top.
        **kwargs:
            This function accepts as argument any of the Index names available associated with a
            single value.
            For example, if 'set_ids' is an available class:`Index <ansys.dpf.post.index.Index>`
            of the class:`DataFrame <ansys.dpf.post.DataFrame>` `df`, then you can plot the data at
            `set_id` 1 by using `df.plot(set_ids=1)`.
            One can get the list of available axes using
            :func:`DataFrame.axes <ansys.dpf.post.DataFrame.axes>`.

        Returns
        -------
            The interactive plotter object used for plotting.

        """
        from ansys.dpf.core.plotter import DpfPlotter as Plotter

        if kwargs != {}:
            self._validate_arguments(arguments=kwargs)
            # Construct the associated label_space
            fc = self.select(**kwargs)._fc
        else:
            # If no kwarg was given, construct a default label_space
            fc = self._fc
        labels = fc.labels
        if "time" in labels:
            label = "time"
            value = fc.get_available_ids_for_label(label)[-1]
            label_space = {label: value}
        elif "frequencies" in labels:
            label = "frequencies"
            value = fc.get_available_ids_for_label(label)[0]
            label_space = {label: value}
        else:
            label_space = fc.get_label_space(0)
        label_space = label_space
        # if "elshape" in self._fc.labels:
        #     merge_solids_shell_op = dpf.operators.logic.solid_shell_fields(fc)
        #     fc = merge_solids_shell_op.eval()

        for field in fc:
            # Treat multi-layer field
            shell_layer_check = field.shell_layers
            if shell_layer_check in [
                shell_layers.topbottom,
                shell_layers.topbottommid,
            ]:
                changeOp = dpf.Operator("change_shellLayers")
                changeOp.inputs.fields_container.connect(fc)
                sl = shell_layers.top
                if shell_layer is not None:
                    if not isinstance(shell_layer, shell_layers):
                        raise TypeError(
                            "shell_layer attribute must be a core.shell_layers instance."
                        )
                    sl = shell_layer
                changeOp.inputs.e_shell_layer.connect(sl.value)  # top layers taken
                fc = changeOp.get_output(0, dpf.types.fields_container)
                break

        fields = fc.get_fields(label_space=label_space)
        plotter = Plotter(**kwargs)
        for field in fields:
            plotter.add_field(field=field, **kwargs)
            # field.plot(text="debug")
        return plotter.show_figure(text=str(label_space), **kwargs)

    def animate(
        self,
        save_as: Union[PathLike, str, None] = None,
        deform: bool = False,
        scale_factor: Union[List[float], float] = 1.0,
        **kwargs,
    ):
        """Animate the result.

        Parameters
        ----------
        save_as : Path of file to save the animation to. Defaults to None. Can be of any format
            supported by pyvista.Plotter.write_frame (.gif, .mp4, ...).
        deform :
            Whether to plot the deformed mesh.
        scale_factor : float, list, optional
            Scale factor to apply when warping the mesh. Defaults to 1.0. Can be a list to make
            scaling frequency-dependent.

        Returns
        -------
            The interactive plotter object used for animation.
        """
        deform_by = None
        if deform:
            try:
                simulation = self._parent_simulation()
                deform_by = simulation._model.results.displacement.on_time_scoping(
                    self._fc.get_time_scoping()
                )
            except Exception as e:
                warnings.warn(
                    UserWarning(
                        "Displacement result unavailable, "
                        f"unable to animate on the deformed mesh:\n{e}"
                    )
                )
        else:
            deform_by = False
        return self._fc.animate(
            save_as=save_as, deform_by=deform_by, scale_factor=scale_factor, **kwargs
        )