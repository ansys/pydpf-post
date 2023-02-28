"""Module containing the ``DataFrame`` class."""
from os import PathLike
from typing import List, Union
import warnings
import weakref

import ansys.dpf.core as dpf

from ansys.dpf.post.index import Index, MultiIndex, location_to_label

display_width = 80
display_max_colwidth = 16


class DataFrame:
    """A DataFrame style API to manipulate DPF data."""

    def __init__(
        self,
        data: Union[dpf.FieldsContainer, None] = None,
        parent_simulation=None,
        index: Union[Index, List[int], None] = None,
        columns: Union[MultiIndex, List[str], None] = None,
    ):
        """Creates a DPF DataFrame based on the given input data.

        Parameters
        ----------
        data:
            Data to use.
        parent_simulation:
            Parent simulation.
        index:
            Index (row labels) to use.
        columns:
            Column labels or class:`ansys.dpf.post.index.MultiIndex` to use.
        """
        self._index = None
        if isinstance(data, dpf.FieldsContainer):
            self._fc = data
            if index is None:
                self._index = Index(
                    name=location_to_label[data[0].location], values=None
                )

        else:
            raise ValueError(
                f"Input data type '{type(data)}' is not a valid data type for DataFrame creation."
            )
        if columns is not None:
            self._columns = columns
        else:
            self._columns = None

        if index is not None:
            self._index = index

        if parent_simulation is not None:
            self._parent_simulation = weakref.ref(parent_simulation)

        self._str = None
        self._last_display_width = display_width
        self._last_display_max_colwidth = display_max_colwidth

        # super().__init__(fields_container._internal_obj, server)

    @property
    def columns(self):
        """Returns the column labels of the DataFrame."""
        if self._columns is None:
            pass
        return self._columns

    @property
    def index(self):
        """Returns the Index for the rows of the DataFrame."""
        if self._index is None:
            pass
        return self._index

    @property
    def axes(self) -> List[str]:
        """Returns a list of the axes of the DataFrame with the row Index and the column Index."""
        names = [self.index.name]
        names.extend(self.columns.label_names)
        names.extend([self.columns.results.name])
        return names

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.FieldsContainer` object."""
        return self._fc

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
        trunc_str = "..."
        # Get the number of rows
        nb_rows = len(self)
        # Get the number of columns
        max_nb_col = width // max_colwidth - 2
        nb_col = len(self.columns)
        if nb_col > max_nb_col:
            max_nb_col = ((width - len(trunc_str)) // max_colwidth - 2) // 2 * 2
            truncate_col = True
        else:
            truncate_col = False

        txt = ""
        # for index_name in self.columns.names:
        #     current_values = getattr(self.columns, index_name).values
        #     mult = nb_col // current_values
        #     txt += index_name.rjust(max_colwidth).join([])

        txt = str(self._fc) + "\n"
        txt += f"DataFrame with columns {self._columns} and index {self._index}"
        txt += "\n"

        self._str = txt

    def plot(self, **kwargs):
        """Plot the result."""
        self._fc[-1].plot(**kwargs)

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
