"""Module containing the ``DataObject`` class."""
from ansys.dpf.post.data import Data
from ansys.dpf.post.errors import PandasImportError


class DataObject:
    """Exposes the fields container generated by a result provider."""

    def __init__(
        self, fields_container=None, server=None, mesh_scoping=None, columns=None
    ):
        """Wrap a FieldsContainer within a DataObject.

        Parameters
        ----------
        fields_container:
            :class:`ansys.dpf.core.fields_container.FieldsContainer`to wrap.
        server:
            DPF server to use.
        mesh_scoping:
            Scoping to use.
        columns:
            Columns to use.
        """
        self._fc = fields_container
        if columns:
            self._columns = columns

        if mesh_scoping:
            self._mesh_scoping = mesh_scoping

        # super().__init__(fields_container._internal_obj, server)

    def __len__(self):
        """Return the length of the DataObject."""
        return len(self._fc)

    def __min__(self, **kwargs):
        """Return the minimum of the data."""
        return self.as_array().min()

    def __max__(self, **kwargs):
        """Return the maximum of the data."""
        return self.as_array().max()

    def __getitem__(self, value):
        """Return a Data in DataObject by index."""
        return Data(self._fc[value])

    def __str__(self):
        """Print DataObject information."""
        txt = f"DPF DataObject:\n {self._fc}"
        return txt

    def max(self, **kwargs):
        """Return the maximum of the data."""
        return float(self.as_array().max())

    def min(self, **kwargs):
        """Return the minimum of the data."""
        return float(self.as_array().min())

    def as_data_frame(self, columns=None, **kwargs):
        """Returns the data from the field container as a Pandas data_frame.

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
        >>> df = displacement.as_data_frame()
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
            columns = self._columns

        # Initialize dataframe
        df = pd.DataFrame()
        for field in self:
            for component in range(field.n_dim):
                column_name = f"{field.name} {columns[component]}"
                data = field.data[:, component]
                ids = field.ids
                df = pd.concat(
                    [df, pd.DataFrame({column_name: data}, index=ids)], axis=1
                )

        # df = pd.DataFrame(
        #     self.as_array(), columns=columns, index=self._mesh_scoping.ids
        # )
        df.name = self._fc[-1].name
        return df

    def as_array(self):
        """Return the DataObject as a NumPy array.

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
        >>> arr = displacement.as_array()
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
