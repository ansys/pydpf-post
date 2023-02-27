"""Module containing the ``Mesh`` class."""
from typing import List

from ansys.dpf.core import MeshedRegion


class Mesh:
    """Exposes the mesh of the solution."""

    def __init__(self, meshed_region: MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    @property
    def available_named_selections(self) -> List[str]:
        """Returns the available named selection of the mesh."""
        return self._meshed_region.available_named_selections

    @property
    def node_ids(self) -> List[int]:
        """Returns the list of node IDs in the mesh."""
        return self._meshed_region.nodes.scoping.ids.tolist()
