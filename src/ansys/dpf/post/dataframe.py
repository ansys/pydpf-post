"""Module containing the ``DataFrame`` class."""
from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar, Union

from ansys.dpf.core import FieldsContainer
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


class DataFrame:
    """A DataFrame style API to manipulate DPF data."""

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
        if columns:
            self.columns = columns
        if index:
            raise NotImplementedError("DataFrame creation using non-DPF entities.")

    def __str__(self) -> str:
        """String representation of the DataFrame."""
        return str(self._fc)

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
