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
        """Returns the available named selection of the mesh."""
        return self._meshed_region.available_named_selections
