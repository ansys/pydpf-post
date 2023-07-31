"""Module containing the ``DataFrame`` class.

DataFrame
---------

"""
from __future__ import annotations

import itertools
from os import PathLike
from typing import List, Union
import warnings

import ansys.dpf.core as dpf
from ansys.dpf.core.dpf_array import DPFArray
from ansys.dpf.core.plotter import DpfPlotter
from ansys.dpf.core.property_fields_container import (
    _MockPropertyFieldsContainer as PropertyFieldsContainer,
)
import ansys.dpf.gate.errors
import numpy as np

from ansys.dpf.post import locations, shell_layers
from ansys.dpf.post.index import (
    CompIndex,
    ElementNodeIndex,
    Index,
    LabelIndex,
    MeshIndex,
    MultiIndex,
    ResultsIndex,
    SetIndex,
    ref_labels,
)

default_display_max_columns = 6
default_display_max_rows = 6


class DataFrame:
    """A DataFrame style API to manipulate DPF data."""

    def __init__(
        self,
        data: Union[dpf.FieldsContainer, PropertyFieldsContainer],
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
        if isinstance(data, dpf.FieldsContainer) or isinstance(
            data, PropertyFieldsContainer
        ):
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

        # If MeshIndex is elemental nodal, update columns
        for i in self.index:
            if isinstance(i, MeshIndex):
                if i.location == locations.elemental_nodal:
                    # Update columns with a dynamic node index
                    self.columns.indexes.append(ElementNodeIndex())

        self._disp_wf = None

        self._str = None
        self._display_max_columns = default_display_max_columns
        self._display_max_rows = default_display_max_rows
        self._last_display_max_columns = self._display_max_columns
        self._last_display_max_rows = self._display_max_rows

        self._last_minmax: dict = {"axis": None, "min": None, "max": None}

    @property
    def display_max_columns(self) -> int:
        """Returns the current maximum number of columns to display for this Dataframe."""
        return self._display_max_columns

    @display_max_columns.setter
    def display_max_columns(self, max_columns: int):
        """Sets a new maximum number of value columns to display for this Dataframe."""
        self._display_max_columns = max_columns

    @property
    def display_max_rows(self) -> int:
        """Returns the current maximum number of rows to display for this Dataframe."""
        return self._display_max_rows

    @display_max_rows.setter
    def display_max_rows(self, max_rows: int):
        """Sets a new maximum number of value rows to display for this Dataframe."""
        self._display_max_rows = max_rows

    @property
    def columns(self) -> MultiIndex:
        """Returns the MultiIndex for the columns of the DataFrame."""
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
    def index(self) -> MultiIndex:
        """Returns the MultiIndex for the rows of the DataFrame."""
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

    @property
    def array(self) -> np.ndarray:
        """Returns the data as a np.ndarray for a single combination of column label values."""
        for index in self.columns:
            if len(index.values) > 1:
                raise ValueError(
                    f"Can only export to array if the DataFrame contains a single {index.name}."
                )
        return self._fc[0].data

    def _filter_arguments(self, arguments):
        """Filter arguments based on available Index names."""
        rows, columns = self.axes
        axes_names = rows.names
        axes_names.extend(columns.names)
        axis_arguments = {}
        keys = list(arguments.keys())
        for argument in keys:
            if argument in axes_names:
                axis_arguments[argument] = arguments.pop(argument)
                # raise ValueError(
                #     f"The DataFrame has no axis {argument}, cannot select it. "
                #     f"Available axes are: {axes_names}."
                # )
        return axis_arguments, arguments

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
        axis_kwargs, _ = self._filter_arguments(arguments=kwargs)
        if ref_labels.set_ids in axis_kwargs.keys():
            axis_kwargs[ref_labels.time] = axis_kwargs[ref_labels.set_ids]
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
        if any([label in axis_kwargs.keys() for label in self._fc.labels]):
            fc = dpf.FieldsContainer()
            to_add_list = []
            for i, field in enumerate(self._fc):
                to_add = True
                label_space = self._fc.get_label_space(i)
                for fc_key in label_space.keys():
                    if fc_key in axis_kwargs.keys() or (
                        fc_key == ref_labels.time
                        and ref_labels.set_ids in axis_kwargs.keys()
                    ):
                        key = fc_key
                        if (
                            fc_key == ref_labels.time
                            and ref_labels.set_ids in axis_kwargs.keys()
                        ):
                            key = ref_labels.set_ids
                        if not isinstance(axis_kwargs[key], list):
                            axis_kwargs[key] = [axis_kwargs[key]]
                        if label_space[fc_key] not in axis_kwargs[key]:
                            to_add = False
                            break
                if to_add:
                    to_add_list.append((label_space, field))
            if len(to_add_list) != 0:
                input_fc = fc
                fc.labels = self._fc.labels
                for label_space, field in to_add_list:
                    fc.add_field(label_space=label_space, field=field)
            else:
                # raise ValueError("Selection criteria resulted in an empty DataFrame.")
                input_fc = fc
                fc.labels = []
                axis_kwargs = {}

        # # Treat selection on components
        if ref_labels.components in axis_kwargs.keys():
            from ansys.dpf.post.simulation import component_label_to_index

            comp_to_extract = axis_kwargs[ref_labels.components]
            if not isinstance(comp_to_extract, list):
                comp_to_extract = [comp_to_extract]
            component_indexes = [component_label_to_index[c] for c in comp_to_extract]
            selector_fc = dpf.operators.logic.component_selector_fc(
                fields_container=input_fc,
                component_number=component_indexes,
                server=server,
            )
            out = selector_fc.outputs.fields_container
            input_fc = out
            comp_index = CompIndex(values=comp_to_extract)

        # Treat selection on DataFrame.index (rescope)
        if isinstance(self.index, MeshIndex):
            mesh_index_name = self.index.name
        else:
            mesh_index_name = self.index.mesh_index.name
        if (
            mesh_index_name in axis_kwargs.keys()
            and mesh_index.location != locations.elemental_nodal
        ):
            if ref_labels.node_ids in mesh_index_name:
                node_ids = axis_kwargs[mesh_index_name]
                if not isinstance(node_ids, (DPFArray, list)):
                    node_ids = [node_ids]
                mesh_scoping = dpf.mesh_scoping_factory.nodal_scoping(
                    node_ids=node_ids,
                    server=server,
                )
            elif ref_labels.element_ids or ref_labels.cell_ids in mesh_index_name:
                if ref_labels.cell_ids in mesh_index_name:
                    location = "cells"
                element_ids = axis_kwargs[mesh_index_name]
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
            if isinstance(input_fc, PropertyFieldsContainer):
                fc = input_fc.rescope(mesh_scoping)
            else:
                rescope_fc = dpf.operators.scoping.rescope_fc(
                    fields_container=input_fc,
                    mesh_scoping=mesh_scoping,
                    server=server,
                )
                out = rescope_fc.outputs.fields_container
                mesh_index = MeshIndex(location=location, values=mesh_scoping.ids)
        elif (
            mesh_index_name in axis_kwargs.keys()
            and mesh_index.location == locations.elemental_nodal
        ):
            raise NotImplementedError(
                "Element selection on a DataFrame with elemental nodal results "
                "is not yet supported"
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

        if "time" in fc.labels:
            set_index = SetIndex(values=fc.get_available_ids_for_label("time"))
        else:
            set_index = SetIndex(values=[])

        column_indexes = [
            results_index,
            set_index,
        ]
        if isinstance(fc, PropertyFieldsContainer):
            column_indexes = [results_index]

        label_indexes = []
        for label in fc.labels:
            if label not in ["time"]:
                fc_values = fc.get_available_ids_for_label(label)
                old_values = self.columns[self.columns.names.index(label)].values
                if old_values is None:
                    values = fc_values
                else:
                    fc_to_old_map = {}
                    for v in old_values:
                        if isinstance(v, str) and ("(" and ")" in v):
                            fc_to_old_map[v.split("(")[1].split(")")[0]] = v
                        else:
                            fc_to_old_map[str(v)] = v
                    values = [fc_to_old_map[str(f)] for f in fc_values]
                column_indexes.append(LabelIndex(name=label, values=values))
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
        axis_kwargs, _ = self._filter_arguments(arguments=kwargs)
        for label in axis_kwargs.keys():
            if label in self.index.names:
                values = getattr(self.index, label).values
            else:
                values = getattr(self.columns, label).values
            indices = axis_kwargs[label]
            if isinstance(indices, list):
                ids = [values[i] for i in indices]
            else:
                ids = values[indices]
            if isinstance(ids, DPFArray):
                ids = ids.tolist()
            axis_kwargs[label] = ids
        return self.select(**axis_kwargs)

    def __len__(self):
        """Return the length of the DataFrame."""
        return len(self.index)

    def __str__(self) -> str:
        """String representation of the DataFrame."""
        if (
            (self._str is None)
            or (self._last_display_max_columns != self.display_max_columns)
            or (self._last_display_max_rows != self.display_max_rows)
        ):
            self._update_str(
                max_columns=self.display_max_columns, max_rows=self.display_max_rows
            )
            self._last_display_max_columns = self.display_max_columns
            self._last_display_max_rows = self.display_max_rows
        return self._str

    def _update_str(self, max_columns: int, max_rows: int):
        """Updates the DataFrame string representation using given display options.

        The string representation is limited to five lines, meaning any DataFrame with more than
        five rows will be truncated and only the five first rows are displayed.

        Parameters
        ----------
        max_columns:
            Maximum number of value columns to show.
        max_rows:
            Maximum number of value rows to show.
        """
        cells = []
        lines = []
        # Static
        empty = ""
        element_node_offset = 0

        # Create cells with row labels and values
        num_column_indexes = len(self.columns)
        num_rows_indexes = len(self.index)

        # Add top-left empty space above row headers
        cells.extend([[empty] * num_column_indexes] * (num_rows_indexes - 1))
        # Add column headers cells
        cells.append([c.name for c in self.columns])
        # Add row headers
        _ = [cells[i].append(row.name) for i, row in enumerate(self._index)]

        # Find out whether rows will be truncated
        truncate_row = 0
        if len(self._fc) > 0:
            num_mesh_entities_to_ask = self._fc[0].size
        else:
            num_mesh_entities_to_ask = 0
        if num_mesh_entities_to_ask > max_rows:
            num_mesh_entities_to_ask = max_rows
            truncate_row = max_rows + num_column_indexes

        comp_values = [1]
        entity_ids = [1]

        # Create row label values combinations
        lists = []
        for index in self.index:
            if isinstance(index, MeshIndex):
                if index._values is not None:
                    values = index.values
                else:
                    values = self._first_n_ids_first_field(num_mesh_entities_to_ask)
                entity_ids = values
            elif isinstance(index, CompIndex):
                values = index.values
                comp_values = values
            else:
                values = index.values
            lists.append(values)
        row_combinations = [p for p in itertools.product(*lists)][:max_rows]

        # Add row labels for the first max_rows combinations
        previous_combination = [None] * len(lists)
        for combination in row_combinations:
            [
                cells[i + element_node_offset].append(
                    str(combination[i])
                    if combination[i] != previous_combination[i]
                    else empty
                )
                for i in range(len(combination))
            ]
            previous_combination = combination

        # Create column label values combinations
        lists = []
        label_positions_in_combinations = {}
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
            if values is None:
                values = [1]
            label_positions_in_combinations[index.name] = position
            lists.append(values)
        column_combinations = [p for p in itertools.product(*lists)]

        # Find out whether columns will be truncated
        truncate_col = 0
        if len(column_combinations) > max_columns:
            truncate_col = max_columns + num_rows_indexes
            column_combinations = column_combinations[:max_columns]

        combination_index = num_rows_indexes
        # For each column combination
        previous_combination = [None] * len(lists)
        for combination in column_combinations[:max_columns]:
            if combination_index > max_columns + num_rows_indexes:
                truncate_col = max_columns + num_rows_indexes
                break
            max_n_col_per_entity = 1
            # ## Fill the label values
            labels = [
                str(combination[i])
                if combination[i] != previous_combination[i]
                else empty
                for i in range(len(combination))
            ]
            previous_combination = combination
            labels.append(empty)  # add empty cell on the row header line
            cells.append(labels)

            # ## Fill the data
            # Create label_space from combination
            label_space = {}
            for label_name in self.labels:
                value = combination[label_positions_in_combinations[label_name]]
                if value is None or value == "":
                    raise ValueError(
                        f"Could not find label value for label {label_name}"
                    )
                elif isinstance(value, str):
                    if "(" in value:
                        value = value.split("(")[1].split(")")[0]
                if label_name == ref_labels.set_ids:
                    label_name = ref_labels.time
                label_space[label_name] = int(value)
            fields = self._fc.get_fields(label_space=label_space)

            # Start counting values found
            n_values = 0
            # Loop over asked entities
            for e, entity_id in enumerate(entity_ids):
                values = []
                # Look for data for this entity in each field
                for field in fields:
                    try:
                        # Try getting data for this entity ID in the current field
                        data = field.get_entity_data_by_id(entity_id).tolist()
                    except ansys.dpf.gate.errors.DPFServerException as e:
                        # If entity data not found in this field, try next field
                        if "Field_GetEntityDataById" in str(e):
                            continue
                        else:
                            raise e
                    except ValueError:
                        continue
                    n_col_per_entity = len(data)
                    # Update max number of node per element in case of elemental nodal
                    if field.location == locations.elemental_nodal:
                        # Update the cells table if more nodes per element than previously
                        if n_col_per_entity > max_n_col_per_entity:
                            for i_n in range(max_n_col_per_entity, n_col_per_entity):
                                to_append = [empty] * (num_column_indexes - 1)
                                to_append.append(str(i_n))
                                to_append.append(empty)
                                to_append.extend(
                                    [empty] * e
                                )  # Pad previous entity lines
                                cells.append(to_append)
                        max_n_col_per_entity = max(
                            max_n_col_per_entity, n_col_per_entity
                        )
                    if n_col_per_entity > max_columns:
                        truncate_col = max_columns + num_rows_indexes
                    # Update number of values found to add per column
                    if isinstance(data[0], list):
                        n_values += len(data[0])
                    else:
                        n_values += 1
                    # # Flatten obtained data and format to string
                    # values_list = [f"{x:.4e}" for y in data for x in y]
                    break_loop = False
                    # Loop over the number of node columns to fill
                    for i_n in range(n_col_per_entity):
                        # Build list of str values to append to current column
                        if i_n < n_col_per_entity:
                            values = data[i_n]
                            if isinstance(values, list):
                                values = [
                                    f"{x}"
                                    if np.issubdtype(type(x), np.integer)
                                    else f"{x:.4e}"
                                    for x in values
                                ]
                            else:
                                values = [
                                    f"{values}"
                                    if np.issubdtype(type(values), np.integer)
                                    else f"{values:.4e}"
                                ]
                        else:
                            values = [empty] * max_rows
                        cells[combination_index + i_n].extend(values)
                        # If already found enough values to print
                        if n_values >= max_rows:
                            # Exit the loop on fields
                            break_loop = True
                    if break_loop:
                        break
                # If already found enough values to print
                if n_values >= max_rows:
                    # Exit the loop on entity IDs
                    break
                if not values:
                    # If no data found for this entity ID, add empty cells
                    values = [[empty] * len(comp_values)]
                    n_values += len(values)
                    cells[combination_index].extend(*values)

            combination_index += max_n_col_per_entity

        self._str = self._format_cells(
            cells=cells, truncate_col=truncate_col, truncate_row=truncate_row
        )

    @staticmethod
    def _format_cells(
        cells: List[List[str]],
        truncate_col: int = 0,
        truncate_row: int = 0,
    ) -> str:
        """Format the data cells into one string.

        Argument cells is a list of column data, each one itself a list of row data.
        """
        truncated_str = "..."

        # Truncate columns
        if truncate_col:
            cells = cells[:truncate_col]
            cells.append([truncated_str] * len(cells[0]))

        # Truncate rows
        if truncate_row:
            for i, _ in enumerate(cells):
                cells[i] = cells[i][: truncate_row + 1]
                cells[i].append(truncated_str)
        n_lines = len(cells[0])

        lines = [""] * n_lines
        for c, column in enumerate(cells):
            if len(column) != n_lines:
                print(cells)
                raise ValueError(
                    f"Column {c} has {len(column)} lines instead of {n_lines}!"
                )
            max_length = max(map(len, column))
            for i in range(n_lines):
                lines[i] += column[i].rjust(max_length + 1)
        lines = [line + "\n" for line in lines]
        return "\n" + "".join(lines)

    def _first_n_ids_first_field(self, n: int) -> List[int]:
        """Return the n first entity IDs of the first field in the underlying FieldsContainer."""
        if len(self._fc) > 0:
            return self._fc[0].scoping.ids[:n]
        else:
            return []

    def __repr__(self):
        """Representation of the DataFrame."""
        return f"DataFrame<index={self.index}, columns={self.columns}>"

    def plot(self, shell_layer=shell_layers.top, **kwargs) -> Union[DpfPlotter, None]:
        """Plot the result.

        Parameters
        ----------
        shell_layer:
            Shell layer to show if multi-layered shell data is present. Defaults to top.
        **kwargs:
            This function accepts as argument any of the Index names available associated with a
            single value.
            For example, if 'set_ids' is an available class:`Index <ansys.dpf.post.index.Index>`
            of the class:`DataFrame <ansys.dpf.post.DataFrame>` `df`, then you can plot the data at
            `set_id` 1 by using `df.plot(set_ids=1)`.
            One can get the list of available axes using
            :func:`DataFrame.axes <ansys.dpf.post.DataFrame.axes>`.
            If the combination of arguments on axes does not return data, this returns None.
            Also supports additional keyword arguments for the plotter. For additional keyword
            arguments, see ``help(pyvista.plot)``.

        Returns
        -------
            The interactive plotter object used for plotting.

        """
        if len(self.index.mesh_index) == 0:
            raise ValueError("Cannot plot a Dataframe with an empty mesh index.")
        label_space = {}
        if kwargs != {}:
            axis_kwargs, kwargs = self._filter_arguments(arguments=kwargs)
            # Construct the associated label_space
            df_temp = self.select(**axis_kwargs)
            fc = df_temp._fc
            if len(fc) == 0:
                warnings.warn(
                    UserWarning(
                        "This combination of criteria did not return data. "
                        f"Nothing to plot for {axis_kwargs}."
                    )
                )
                return
            labels = fc.labels
            # TODO: treat complex label by taking amplitude
            for label in labels:
                if label == "time":
                    value = fc.get_available_ids_for_label(label)[-1]
                elif label == "frequencies":
                    value = fc.get_available_ids_for_label(label)[0]
                elif label in axis_kwargs.keys():
                    value = axis_kwargs[label]
                    if isinstance(value, list):
                        raise ValueError(
                            f"Plot argument '{label}' must be a single value."
                        )
                else:
                    value = fc.get_available_ids_for_label(label)[0]
                label_space[label] = value
        else:
            axis_kwargs = {}
            # If no kwarg was given, construct a default label_space
            fc = self._fc
            label_space = fc.get_label_space(0)

        if len(fc) == 0:
            raise ValueError("No data to plot.")
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
                            "shell_layer attribute must be a dpf.shell_layers instance."
                        )
                    sl = shell_layer
                changeOp.inputs.e_shell_layer.connect(sl.value)  # top layers taken
                fc = changeOp.get_output(0, dpf.types.fields_container)
                break

        # Merge fields for all 'elshape' label values if none selected
        if "elshape" in self._fc.labels and "elshape" not in axis_kwargs.keys():
            merge_solids_shell_op = dpf.operators.logic.solid_shell_fields(fc)
            fc = merge_solids_shell_op.eval()

        # Merge fields for all 'stage' label values if none selected
        if "stage" in self._fc.labels and "stage" not in axis_kwargs.keys():
            merge_stages_op = dpf.operators.utility.merge_fields_by_label(
                fields_container=fc, label="stage"
            )
            fc = merge_stages_op.outputs.fields_container()
            label_space.pop("stage")

        fields = fc.get_fields(label_space=label_space)
        # for field in fields:
        if len(fields) > 1:
            # try:
            #     for field in fields[1:]:
            #         plotter.add_field(field=field, **kwargs)
            # except Exception as e:
            raise ValueError(
                f"Plotting failed with filter {axis_kwargs} due to incompatible data."
            )
            # warnings.warn(
            #     UserWarning(
            #         "Plotting criteria resulted in incompatible data. "
            #         "Try narrowing down to specific values for each column."
            #     )
            # )
            # return None
        # for field in fields:
        if len(fields) > 1:
            # try:
            #     for field in fields[1:]:
            #         plotter.add_field(field=field, **kwargs)
            # except Exception as e:
            raise ValueError(
                f"Plotting failed with filter {axis_kwargs} due to incompatible data."
            )
            # warnings.warn(
            #     UserWarning(
            #         "Plotting criteria resulted in incompatible data. "
            #         "Try narrowing down to specific values for each column."
            #     )
            # )
            # return None
        # field.plot(text="debug")
        field_to_plot = fields[0]
        # If multi-component, take the norm
        if field_to_plot.component_count > 1:
            field_to_plot = dpf.operators.math.norm(
                field_to_plot, server=field_to_plot._server
            ).eval()
        plotter = DpfPlotter(**kwargs)
        plotter.add_field(field=field_to_plot, **kwargs)

        return plotter.show_figure(
            title=kwargs.pop("title", str(label_space)), **kwargs
        )

    def animate(
        self,
        save_as: Union[PathLike, str, None] = None,
        deform: bool = False,
        scale_factor: Union[List[float], float] = 1.0,
        shell_layer: shell_layers = shell_layers.top,
        **kwargs,
    ):
        """Animate the DataFrame along its 'set' axis.

        .. note::
            At the moment only useful to produce a temporal animation. Each frame will correspond
            to data for a value of the SetIndex in the columns MultiIndex.

        Parameters
        ----------
        save_as : Path of file to save the animation to. Defaults to None. Can be of any format
            supported by pyvista.Plotter.write_frame (.gif, .mp4, ...).
        deform :
            Whether to plot the deformed mesh.
        scale_factor : float, list, optional
            Scale factor to apply when warping the mesh. Defaults to 1.0. Can be a list to make
            scaling frequency-dependent.
        shell_layer:
            Shell layer to show if multi-layered shell data is present. Defaults to top.

        Returns
        -------
            The interactive plotter object used for animation.
        """
        deform_by = None
        if deform:
            if self._disp_wf is None:
                warnings.warn(
                    UserWarning(
                        "Displacement result unavailable, "
                        "unable to animate on the deformed mesh."
                    )
                )
            else:
                wf = dpf.workflow.Workflow()
                forward_op = dpf.operators.utility.forward_fields_container(
                    server=self._fc._server
                )
                wf.add_operator(forward_op)
                wf.set_input_name("input", forward_op.inputs.fields)
                output_input_names = ("output", "input")
                wf.connect_with(
                    left_workflow=self._disp_wf, output_input_names=output_input_names
                )

                deform_by = forward_op
        else:
            deform_by = False

        fc = self._fc

        # Modify fc to merge fields of different eltype at same time and to select a shell_layer
        # (until it is done on core-side -> animation feature to refactor)

        sl = shell_layers.top
        # Select shell_layer
        for field in fc:
            # Treat multi-layer field
            shell_layer_check = field.shell_layers
            if shell_layer_check in [
                shell_layers.topbottom,
                shell_layers.topbottommid,
            ]:
                if shell_layer is not None:
                    if not isinstance(shell_layer, shell_layers):
                        raise TypeError(
                            "shell_layer attribute must be a dpf.shell_layers instance."
                        )
                    sl = shell_layer
                shell_layer_op = dpf.operators.utility.change_shell_layers(
                    fields_container=fc,
                    server=fc._server,
                    e_shell_layer=sl.value(),
                )
                fc = shell_layer_op.outputs.fields_container_as_fields_container()
                break

        # Set shell layer input for self._disp_wf
        self._disp_wf.connect("shell_layer_int", sl.value)

        # Merge shell and solid fields at same set
        if "eltype" in fc.labels:
            merge_op = dpf.operators.utility.merge_fields_by_label(
                fields_container=fc,
                label="eltype",
            )
            fc = merge_op.outputs.fields_container()

        return fc.animate(
            save_as=save_as, deform_by=deform_by, scale_factor=scale_factor, **kwargs
        )

    def min(self, axis: Union[int, str, None] = 0) -> Union[DataFrame, float]:
        """Return the minimum value over the requested axis.

        Parameters
        ----------
        axis:
          Axis to perform minimum across.
          Defaults to the MeshIndex (0), the row index containing mesh entity IDs.
          This computes the minimum across the mesh for each set.
          Can also be the SetIndex (1), the column index containing set (time/frequency) IDs.
          This computes the minimum across sets (time/frequency) for each mesh entity.

        Returns
        -------
          A scalar if the result of the query is a single number,
          or a DataFrame if several numbers along one or several axes.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.StaticMechanicalSimulation(examples.download_crankshaft())
        >>> displacement = simulation.displacement(all_sets=True)
        >>> # Compute the maximum displacement value for each component at each time-step
        >>> minimum_over_mesh = displacement.min(axis="node_ids")
        >>> print(minimum_over_mesh)  # doctest: +NORMALIZE_WHITESPACE
           results       U (m)
           set_ids           1           2           3
        components
                 X -7.4732e-04 -1.5081e-03 -2.2755e-03
                 Y -4.0138e-04 -8.0316e-04 -1.2014e-03
                 Z -2.1555e-04 -4.3299e-04 -6.5101e-04
        >>> # Compute the maximum displacement for each node and component across time
        >>> minimum_over_time = displacement.min(axis="set_ids")
        >>> print(minimum_over_time)  # doctest: +NORMALIZE_WHITESPACE
                    results       U (m)
         node_ids components
             4872          X -3.4137e-05
                           Y  5.1667e-04
                           Z -4.1346e-06
             9005          X -5.5625e-05
                           Y  4.8445e-04
                           Z -4.9795e-07
              ...        ...         ...
        >>> # Compute the maximum displacement overall
        >>> minimum_overall = minimum_over_time.min()
        >>> print(minimum_overall)  # doctest: +NORMALIZE_WHITESPACE
           results       U (m)
        components
                 X -2.2755e-03
                 Y -1.2014e-03
                 Z -6.5101e-04
        """
        self._query_min_max(axis)
        return self._last_minmax["min"]

    def max(self, axis: Union[int, str, None] = 0) -> Union[DataFrame, float]:
        """Return the maximum value over the requested axis.

        Parameters
        ----------
        axis:
          Axis to perform maximum across.
          Defaults to the MeshIndex (0), the row index containing mesh entity IDs.
          This computes the maximum across the mesh for each set.
          Can also be the SetIndex (1), the column index containing set (time/frequency) IDs.
          This computes the maximum across sets (time/frequency) for each mesh entity.

        Returns
        -------
          A scalar if the result of the query is a single number,
          or a DataFrame if several numbers along one or several axes.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.StaticMechanicalSimulation(examples.download_crankshaft())
        >>> displacement = simulation.displacement(all_sets=True)
        >>> # Compute the maximum displacement value for each component at each time-step
        >>> maximum_over_mesh = displacement.max(axis="node_ids")
        >>> print(maximum_over_mesh)  # doctest: +NORMALIZE_WHITESPACE
           results       U (m)
           set_ids           1           2           3
        components
                 X  7.3303e-04  1.4495e-03  2.1441e-03
                 Y  1.3962e-03  2.7884e-03  4.1656e-03
                 Z  2.1567e-04  4.3321e-04  6.5135e-04
        >>> # Compute the maximum displacement for each node and component across time
        >>> maximum_over_time = displacement.max(axis="set_ids")
        >>> print(maximum_over_time)  # doctest: +NORMALIZE_WHITESPACE
                     results       U (m)
         node_ids components
             4872          X  5.6781e-06
                           Y  1.5417e-03
                           Z -2.6398e-06
             9005          X -2.6323e-06
                           Y  1.4448e-03
                           Z  5.3134e-06
              ...        ...         ...
        >>> # Compute the maximum displacement overall
        >>> maximum_overall = maximum_over_time.max()
        >>> print(maximum_overall)  # doctest: +NORMALIZE_WHITESPACE
           results       U (m)
        components
                 X  2.1441e-03
                 Y  4.1656e-03
                 Z  6.5135e-04
        """
        self._query_min_max(axis)
        return self._last_minmax["max"]

    def _query_min_max(self, axis: Union[int, str, None]) -> None:
        """Create a DPF workflow based on the query arguments for min/max."""
        # Translate None query to empty dict
        if axis in [None, 0, self.index.mesh_index.name]:
            axis = 0
        elif axis in [1, ref_labels.set_ids]:
            axis = 1
        else:
            raise ValueError(f"'{axis}' is not an available axis value.")

        # print(f"{axis=}")
        # If same query as last and last is not None, do not change
        if self._last_minmax["axis"] == axis and not self._last_minmax["axis"] is None:
            return
        # If in need of an update, create the appropriate workflow
        wf = dpf.Workflow(server=self._fc._server)
        wf.progress_bar = False

        # If over mesh
        if axis == 0:
            min_max_op = dpf.operators.min_max.min_max_over_label_fc(
                fields_container=self._fc,
                label="time",
                server=self._fc._server,
            )
            # Here the fields are located on the label ("time"), so we have to "transpose" it.
            # Extract the data for each time (entity) from the field and create a fields_container

            min_fc = dpf.FieldsContainer(server=self._fc._server)
            min_fc.add_label(label="time")
            min_field = min_max_op.outputs.field_min()
            for i, time in enumerate(min_field.scoping.ids):
                min_fc.add_field(
                    label_space={"time": time},
                    field=dpf.fields_factory.field_from_array(
                        arr=min_field.get_entity_data(i),
                        server=self._fc._server,
                    ),
                )

            max_fc = dpf.FieldsContainer(server=self._fc._server)
            max_fc.add_label(label="time")
            max_field = min_max_op.outputs.field_max()
            for i, time in enumerate(max_field.scoping.ids):
                max_fc.add_field(
                    label_space={"time": time},
                    field=dpf.fields_factory.field_from_array(
                        arr=max_field.get_entity_data(i),
                        server=self._fc._server,
                    ),
                )

            index = MultiIndex(
                indexes=[i for i in self.index if i != self.index.mesh_index]
            )
            columns = self.columns

        # If over time
        else:
            min_max_op = dpf.operators.min_max.min_max_over_time_by_entity(
                fields_container=self._fc,
                server=self._fc._server,
            )
            wf.set_output_name("min", min_max_op.outputs.min)
            wf.set_output_name("max", min_max_op.outputs.max)

            index = self.index
            columns = MultiIndex(
                indexes=[c for c in self.columns if c != self.columns.set_ids]
            )

            min_fc = wf.get_output("min", dpf.types.fields_container)
            max_fc = wf.get_output("max", dpf.types.fields_container)

        self._last_minmax["min"] = DataFrame(
            data=min_fc,
            index=index,
            columns=columns,
        )
        self._last_minmax["max"] = DataFrame(
            data=max_fc,
            index=index,
            columns=columns,
        )
        self._last_minmax["axis"] = axis
