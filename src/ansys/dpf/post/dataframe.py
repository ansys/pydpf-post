"""Module containing the ``DataFrame`` class."""
from __future__ import annotations

from typing import TYPE_CHECKING, List, TypeVar, Union

from ansys.dpf.core import FieldsContainer, ScopingsContainer
from ansys.dpf.core.operators.utility import merge_scopings
import numpy
import numpy as np

from ansys.dpf.post.errors import PandasImportError

PANDAS_PRESENT = False
if TYPE_CHECKING:
    try:
        import pandas as pd

        PANDAS_PRESENT = True
    except ImportError:
        pass

if PANDAS_PRESENT:
    PandasDataFrameType = pd.DataFrame
else:
    PandasDataFrameType = TypeVar("PandasDataFrameType")

display_width = 80
display_max_colwidth = 16


class DataFrame:
    """A DataFrame style API to manipulate DPF data."""

    trunc_str = "..."

    def __init__(
        self,
        data: Union[FieldsContainer, numpy.ndarray, PandasDataFrameType],
        index=None,
        columns=None,
        server=None,
    ):
        """Creates a DPF DataFrame based on the given input data.

        Parameters
        ----------
        data:
            Data to initialize the dataframe with.
        index:
            Sequencing and Mesh multi-index to use when instantiating from a numpy.ndarray
            or a pandas.DataFrame.
        columns:
            Columns labels to use.
        server:
            DPF server to use.
        """
        self._server = server
        # When created from a DPF FieldsContainer
        if type(data) == FieldsContainer:
            # Simply store the FieldsContainer
            self._fc = data
        # When created from non-DPF data structures
        else:
            # We need to instantiate a corresponding FieldsContainer server-side
            if type(data) == np.ndarray:
                raise NotImplementedError("DataFrame creation using np.ndarray.")
            elif type(data) == PandasDataFrameType:
                raise NotImplementedError("DataFrame creation using pd.DataFrame.")
            else:
                raise ValueError(
                    "Argument 'data' must be either a numpy.ndarray, "
                    "a pandas.DataFrame or an ansys.dpf.core.FieldsContainer."
                )
        self._columns = columns
        if index:
            raise NotImplementedError("DataFrame creation using non-DPF entities.")

        self._str = None
        self._len = None
        self._last_display_width = display_width
        self._last_display_max_colwidth = display_max_colwidth

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

    def head(self, n: int = 5) -> str:
        """String representation of the n first rows of the DataFrame."""
        nb_rows = len(self)
        if nb_rows < n:
            n = nb_rows
        return self._print_rows(indexes=list(range(0, n)))

    def tail(self, n: int = 5) -> str:
        """String representation of the n last rows of the DataFrame."""
        nb_rows = len(self)
        if nb_rows < n:
            n = nb_rows
        return self._print_rows(indexes=list(range(nb_rows - n, nb_rows)))

    def _print_rows(
        self,
        indexes: List,
    ) -> str:
        rows_str, truncate_col, col_indexes, col_indexes_next = self._build_columns_str(
            display_width, display_max_colwidth
        )
        for i in indexes:
            rows_str += self._build_row_str(
                i, display_max_colwidth, truncate_col, col_indexes, col_indexes_next
            )
        return rows_str

    def _build_row_str(
        self,
        row_index: int,
        max_colwidth: int,
        truncate_col: bool,
        col_indexes: List[int],
        col_indexes_next: List[int],
    ):
        set_id, mesh_id, values = self._get_row_data(row_index)
        row_str = str(set_id).rjust(max_colwidth)
        row_str += str(mesh_id).rjust(max_colwidth)
        if len(values) > 0:
            if type(values[0]) not in [float, int]:
                values = [v for value in values for v in value]
            row_str += "".join(
                [f"{values[i]:.6e}".rjust(max_colwidth) for i in col_indexes]
            )
            if truncate_col:
                row_str += "".join(self.trunc_str).join(
                    [f"{values[i]:.6e}".rjust(max_colwidth) for i in col_indexes_next]
                )
        return row_str + "\n"

    def _build_columns_str(self, width, max_colwidth):
        # Get the number of columns
        max_nb_col = width // max_colwidth - 2
        # print("max_nb_col", max_nb_col)
        nb_col = len(self.columns)
        if nb_col > max_nb_col:
            max_nb_col = ((width - len(self.trunc_str)) // max_colwidth - 2) // 2 * 2
            col_indexes = list(range(max_nb_col // 2))
            col_indexes_next = list(range(nb_col - (max_nb_col // 2), nb_col))
            truncate_col = True
        else:
            col_indexes = list(range(nb_col))
            truncate_col = False
            col_indexes_next = []
        labels = ["set ID", "node ID"]
        labels.extend([self.columns[i] for i in col_indexes])
        if truncate_col:
            labels.extend(self.trunc_str)
            labels.extend([self.columns[i] for i in col_indexes_next])
        return (
            "".join([label.rjust(max_colwidth) for label in labels]) + "\n",
            truncate_col,
            col_indexes,
            col_indexes_next,
        )

    def _update_str(self, width: int, max_colwidth: int):
        """Updates the DataFrame string representation using given display options.

        The string representation is limited to five lines, meaning any DataFrame with more than
        five rows will be truncated and only the two first and two last rows are displayed.

        Parameters
        ----------
        width:
            Number of characters to use for the total width.
        max_colwidth:
            Maximum number of characters to use for each column.
        """
        # Get the number of rows
        nb_rows = len(self)
        # print(f"{nb_rows=}")
        if nb_rows > 5:
            truncate_rows = True
            row_indexes = [0, 1]
            row_indexes_next = [nb_rows - 2, nb_rows - 1]
        else:
            truncate_rows = False
            row_indexes = list(range(nb_rows))

        txt, truncate_col, col_indexes, col_indexes_next = self._build_columns_str(
            width, max_colwidth
        )

        for row_index in row_indexes:
            txt += self._build_row_str(
                row_index, max_colwidth, truncate_col, col_indexes, col_indexes_next
            )

        if truncate_rows:
            txt += "".join(
                [
                    self.trunc_str.rjust(max_colwidth)
                    for _ in range(len(col_indexes) + 2)
                ]
            )
            if truncate_col:
                txt += self.trunc_str.rjust(max_colwidth).join(
                    [self.trunc_str.rjust(max_colwidth) for _ in col_indexes_next]
                )
            txt += "\n"
            for row_index in row_indexes_next:
                txt += self._build_row_str(
                    row_index, max_colwidth, truncate_col, col_indexes, col_indexes_next
                )

        self._str = txt

    def _get_row_data(self, row_index):
        # Find the set and entity ID corresponding to this row
        set_index_in_fc = row_index // self._n_unique_mesh_entities
        mesh_index = row_index % len(self._unique_mesh_ids)
        mesh_id = self._unique_mesh_ids[mesh_index]
        time_id = self._fc.get_label_scoping()[set_index_in_fc]
        # Get the data
        fields = self._fc.get_fields_by_time_complex_ids(timeid=time_id)
        # Build the str representation
        values = []
        for field in fields:
            try:
                values.extend(field.get_entity_data_by_id(mesh_id))
                # print(f"{values=}")
            except Exception as e:
                raise e
        return time_id, mesh_id, values

    def __len__(self):
        """Returns the number of unique data sets in the DataFrame."""
        if self._len is None:
            self._update_len()
        return self._len

    def _update_len(self):
        """Compute and store the DataFrame row length from underlying Fields."""
        scoping_list = [f.scoping for f in self._fc]
        sc = ScopingsContainer(server=self._fc._server)
        sc.labels = ["field"]
        for i, scoping in enumerate(scoping_list):
            sc.add_scoping({"field": i}, scoping)
        global_scoping = merge_scopings(scopings=sc, server=self._fc._server).eval()
        n_unique_mesh_ids = global_scoping.size
        self._unique_mesh_ids = global_scoping.ids
        self._len = n_unique_mesh_ids

    @property
    def columns(self) -> list:
        """Returns a list of column labels."""
        if self._columns is None:
            self._update_columns()
        return self._columns

    def _update_columns(self):
        """Update the list of columns from the underlying FieldsContainer."""
        columns = []
        for field in self._fc:
            result = field.name
            nb_components = field.component_count
            for comp in range(nb_components):
                columns.append(result.split("_")[0] + str(comp))
        self._columns = list(set(columns))

    def to_pandas(self, columns=None, **kwargs) -> PandasDataFrameType:
        """Returns the current DPF DataFrame as a Pandas DataFrame.

        Parameters
        ----------
        columns :
            Columns description.

        Returns
        -------
        class:`pandas.core.frame.DataFrame`

        Examples
        --------
        >>> import pandas as pd
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.multishells_rst)
        >>> # Export the displacements vector field at step 1 as a DataFrame
        >>> displacement = simulation.displacement(steps=[1], nodes=[1, 2, 3])
        >>> df = displacement.to_pandas()
        >>> print(df)
                  X          Y         Z
        1  0.398320 -13.797378 -0.163767
        2  0.233114 -13.797652 -0.153190
        3  0.367542 -13.808151 -0.163739
        >>> print(df.name)
        displacement_1.s

        """
        try:
            import pandas as pd
        except ModuleNotFoundError:
            raise PandasImportError

        if not columns:
            columns = self.columns
        # The wrapped FieldsContainer should only contain one Field
        if len(self) > 1:
            raise NotImplementedError(
                "This function requires for the DataFrame to contain only"
                "one data entity."
            )
        df = pd.DataFrame(self.to_numpy(), columns=columns, index=self.index.mesh.ids)
        df.name = self._fc[-1].name
        return df

    def to_numpy(self):
        """Return the current DPF DataFrame as a NumPy ndarray.

        Returns
        -------
        class:`numpy.ndarray`

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.multishells_rst)
        >>> # Export the displacements vector field at step 1 as a numpy.ndarray
        >>> displacement = simulation.displacement(steps=[1], nodes=[1, 2, 3])
        >>> arr = displacement.to_numpy()
        >>> print(arr)
        [[  0.39831985 -13.79737819  -0.16376683]
         [  0.23311378 -13.79765179  -0.15318972]
         [  0.36754181 -13.80815135  -0.16373945]]

        """
        return self._fc[-1].data

    def plot(self, **kwargs):
        """Plot the result."""
        self._fc[-1].plot(**kwargs)

    def animate(self, **kwargs):
        """Animate the result.

        Returns
        -------
            The interactive plotter object used for animation.
        """
        return self._fc.animate(**kwargs)
