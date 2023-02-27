"""Module containing the ``DataFrame`` class."""
from os import PathLike
from typing import List, Union
import warnings
import weakref


class DataFrame:
    """A DataFrame style API to manipulate DPF data."""

    def __init__(
        self,
        fields_container=None,
        parent_simulation=None,
        index=None,
        columns=None,
    ):
        """Creates a DPF DataFrame based on the given input data.

        Parameters
        ----------
        fields_container:
            :class:`ansys.dpf.core.fields_container.FieldsContainer`to wrap.
        parent_simulation:
            Parent simulation.
        index:
            Index to use.
        columns:
            Columns labels to use.
        """
        self._fc = fields_container
        if columns is not None:
            self._columns = columns

        if index is not None:
            self._index = index

        if parent_simulation is not None:
            self._parent_simulation = weakref.ref(parent_simulation)

        # super().__init__(fields_container._internal_obj, server)

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.FieldsContainer` object."""
        return self._fc

    def __len__(self):
        """Return the length of the DataFrame."""
        return len(self._fc)

    def __str__(self) -> str:
        """String representation of the DataFrame."""
        return str(self._fc)

    def __min__(self, **kwargs):
        """Return the minimum of the data."""
        return self.as_array().min()

    def __max__(self, **kwargs):
        """Return the maximum of the data."""
        return self.as_array().max()

    def max(self, **kwargs):
        """Return the maximum of the data."""
        return float(self.as_array().max())

    def min(self, **kwargs):
        """Return the minimum of the data."""
        return float(self.as_array().min())

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
                pass
        return self._fc.animate(
            save_as=save_as, deform_by=deform_by, scale_factor=scale_factor, **kwargs
        )
