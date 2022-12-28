"""Module containing the ``Mesh`` class."""
from typing import List

from ansys.dpf.core import MeshedRegion


class Mesh:
    """Exposes the mesh of the simulation.

    Parameters
    ----------
    meshed_region:
        Instance of :class:`MeshedRegion <ansys.dpf.core.meshed_region.MeshedRegion>` to wrap.

    """

    def __init__(self, meshed_region: MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    @property
    def available_named_selections(self) -> List[str]:
        """Return the available named selection of the mesh."""
        return self._meshed_region.available_named_selections

    @property
    def available_property_fields(self) -> List[str]:
        """Return the available property fields of the mesh."""
        return self._meshed_region.available_property_fields

    @property
    def grid(self):
        """Return the grid of the mesh."""
        return self._meshed_region.grid

    @property
    def nodes(self):
        """Return the nodes of the mesh."""
        return self._meshed_region.nodes

    @property
    def elements(self):
        """Return the elements of the mesh."""
        return self._meshed_region.elements

    @property
    def unit(self):
        """Return the unit of the mesh."""
        return self._meshed_region.unit
