"""Module containing the ``Data`` class."""


class Data:
    """Wraps the PyDPF-Core :class:`ansys.dpf.core.field.Field` object.

    Parameters
    ----------
    TODO

    Examples
    --------
    TODO

    """

    def __init__(self, field):
        """Initialize Data class."""
        self._field = field

    @property
    def name(self):
        """Return name of Data."""
        return self._field.name

    @property
    def location(self):
        """Return location of Data."""
        return self._field.location

    @property
    def data(self):
        """Return data array of Data."""
        return self._field.data

    @property
    def ids(self):
        """Return IDs of Data."""
        return self._field.scoping.ids

    @property
    def n_dim(self):
        """Return number of dimensions of Data."""
        return self._field.component_count

    @property
    def n_data(self):
        """Return number of entries in Data.data."""
        return self._field.elementary_data_count

    @property
    def time_freq(self):
        """Return time/frequencies of Data."""
        return self._field.time_freq_support.time_frequencies

    @property
    def unit(self):
        """Return unit of Data."""
        return self._field.unit

    def __str__(self):
        """Print Data information."""
        from ansys.dpf.core.core import _description

        return _description(self._field._internal_obj, self._field._server)
