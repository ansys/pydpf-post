"""Module containing the ``Mesh`` class."""

from ansys.dpf.core import MeshedRegion


class Mesh:
    """Exposes the mesh of the solution."""

    def __init__(self, meshed_region: MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    def __str__(self):
        """String representation of this class."""
        return str(self._meshed_region)

    @property
    def available_named_selections(self):
        """Returns the available named selection of the mesh."""
        return self._meshed_region.available_named_selections
