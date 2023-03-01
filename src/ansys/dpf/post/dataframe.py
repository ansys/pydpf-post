"""Module containing the ``DataFrame`` class."""
from __future__ import annotations

import itertools
from os import PathLike
from typing import List, Union
import warnings

import ansys.dpf.core as dpf
from ansys.dpf.core.dpf_array import DPFArray

from ansys.dpf.post.index import (
    CompIndex,
    FrequencyIndex,
    Index,
    LabelIndex,
    MeshIndex,
    MultiIndex,
    ResultsIndex,
    SetIndex,
    TimeIndex,
)

display_width = 80
display_max_colwidth = 10


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
    def axes(self) -> List[str]:
        """Returns a list of the axes of the DataFrame with the row Index and the column Index."""
        names = self.index.names
        names.extend(self.columns.names)
        return names

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
        """Returns a list of the names of available ResultsIndex indexes."""
        names = [index.name for index in self.index if isinstance(index, LabelIndex)]
        names.extend(
            [index.name for index in self.columns if isinstance(index, LabelIndex)]
        )
        return names

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.FieldsContainer` object."""
        return self._fc

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
        # Check for invalid arguments
        axes = self.axes
        for argument in kwargs.keys():
            if argument not in axes:
                raise ValueError(
                    f"The DataFrame has no axis {argument}, cannot select it. "
                    f"Available axes are: {axes}."
                )
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
                for key in label_space.keys():
                    if not isinstance(kwargs[key], list):
                        kwargs[key] = [kwargs[key]]
                    if label_space[key] in kwargs[key]:
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
        if mesh_index_name in kwargs.keys():
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
            `time` value by using `df.select(time=0)`.
            One can get the list of available axes using
            :func:`DataFrame.axes <ansys.dpf.post.DataFrame.axes>`.

        Returns
        -------
            A DataFrame of the selected values.

        """
        # Check for invalid arguments
        axes = self.axes
        for argument in kwargs.keys():
            if argument not in axes:
                raise ValueError(
                    f"The DataFrame has no axis {argument}, cannot select it. "
                    f"Available axes are: {axes}."
                )
        for label in kwargs.keys():
            indices = kwargs[label]
            if label in self.index.names:
                ids = getattr(self.index, label).values[indices]
            else:
                ids = getattr(self.columns, label).values[indices]
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
        empty = " " * max_colwidth
        truncated = "...".rjust(max_colwidth)
        max_n_col = width // max_colwidth
        max_n_rows = 5
        # Create lines with row labels and values
        num_column_indexes = len(self.columns)
        num_rows_indexes = len(self.index)
        n_max_value_col = max_n_col - num_rows_indexes
        column_headers = []
        for column_index in self.columns:
            column_headers.append(
                empty * (num_rows_indexes - 1) + column_index.name.rjust(max_colwidth)
            )

        lines = column_headers

        row_headers = "".join(
            [row_index.name.rjust(max_colwidth) for row_index in self.index]
        )
        lines.append(row_headers)
        entity_ids = []
        # Create combinations for rows
        num_mesh_entities_to_ask = 5
        lists = []
        entity_ids = None
        for index in self.index:
            if isinstance(index, MeshIndex):
                if index._values is not None:
                    values = index.values
                    entity_ids = values
                else:
                    values = self._first_n_ids_first_field(num_mesh_entities_to_ask)
            elif isinstance(index, CompIndex):
                values = index.values
            else:
                values = index.values
            lists.append(values)
        row_combinations = [p for p in itertools.product(*lists)]

        # Add row headers for the first combinations (max_n_rows)
        previous_combination = [None] * len(lists)
        for combination in row_combinations[:max_n_rows]:
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
        num_mesh_entities_to_ask = 5
        entity_ids = None
        lists = []
        time_position = 1
        complex_position = None
        for position, index in enumerate(self.columns):
            if isinstance(index, MeshIndex):
                if index._values is not None:
                    values = index.values
                    entity_ids = values
                else:
                    values = self._first_n_ids_first_field(num_mesh_entities_to_ask)
            elif isinstance(index, CompIndex):
                values = index.values
            elif isinstance(index, (FrequencyIndex, TimeIndex)):
                values = index.values
                time_position = position
            else:
                values = index.values
                if index.name == "complex":
                    complex_position = position
            lists.append(values)
        combinations = [p for p in itertools.product(*lists)]

        # Query data on entity_ids

        # Add text for the first n_max_value_col columns
        previous_combination = [None] * len(lists)
        for combination in combinations[:n_max_value_col]:
            to_append = [
                str(combination[i]).rjust(max_colwidth)
                if combination[i] != previous_combination[i]
                else empty
                for i in range(len(combination))
            ]
            to_append.append(empty)
            # Get data in the FieldsContainer for those positions
            # Create label_space from combination
            label_space = {"time": combination[time_position]}
            if complex_position is not None:
                label_space["complex"] = combination[complex_position]
            fields = self._fc.get_fields(label_space=label_space)
            values = []
            if entity_ids is None:
                for field in fields:
                    array_values = []
                    for k in list(range(num_mesh_entities_to_ask)):
                        try:
                            values_list = field.get_entity_data(k).tolist()
                            array_values.append(values_list)
                        except Exception as e:
                            pass
                    array_values = [
                        item for sublist in array_values for item in sublist
                    ]
                    array_values = [
                        item for sublist in array_values for item in sublist
                    ]
                    values.extend(array_values)
            else:
                for entity_id in entity_ids:
                    for field in fields:
                        array_values = field.get_entity_data_by_id(entity_id).tolist()
                        if array_values != []:
                            array_values = [
                                item for sublist in array_values for item in sublist
                            ]
                            array_values = [
                                item for sublist in array_values for item in sublist
                            ]
                            values.extend(array_values)

            value_strings = [
                f"{value:.{max_colwidth - 8}e}".rjust(max_colwidth)
                for value in values[: len(lines) - (num_column_indexes + 1)]
                # TODO: error here, must choose the correct components
            ]
            to_append.extend(value_strings)
            previous_combination = combination
            # print(to_append)
            # print(len(to_append))
            # print(len(lines))
            for i in range(len(lines)):
                lines[i] = lines[i] + to_append[i]

        txt = "\n" + "".join([line + "\n" for line in lines])
        self._str = txt

    def _first_n_ids_first_field(self, n: int):
        """Return the n first entity IDs of the first field in the underlying FieldsContainer."""
        return self._fc[0].scoping.ids[:n]

    def __repr__(self):
        """Representation of the DataFrame."""
        return f"DataFrame<index={self.index}, columns={self.columns}>"

    def plot(self, **kwargs):
        """Plot the result.

        Parameters
        ----------
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
        if kwargs != {}:
            # Check for invalid arguments
            axes = self.axes
            for argument in kwargs.keys():
                if argument not in axes:
                    raise ValueError(
                        f"The DataFrame has no axis {argument}, cannot plot it. "
                        f"Available axes are: {axes}."
                    )
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
        field = fc.get_field(label_space_or_index=label_space)
        return field.plot(text=str(label_space), **kwargs)

    def animate(
        self,
        save_as: Union[PathLike, None] = None,
        deform: bool = False,
        scale_factor: Union[List[float], float, None] = None,
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
        return self._fc.animate(
            save_as=save_as, deform_by=deform_by, scale_factor=scale_factor, **kwargs
        )