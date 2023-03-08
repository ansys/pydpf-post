"""Module containing the ``Mesh`` class.

Mesh
----

"""
from typing import List

from ansys.dpf.core import MeshedRegion


class Mesh:
    """Exposes the mesh of the solution."""

    def __init__(self, meshed_region: MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    def __str__(self):
        """String representation of this class."""
        return str(self._meshed_region).replace("Meshed Region", "Mesh")

    @property
    def available_named_selections(self) -> List[str]:
        """Returns the available named selections of the mesh."""
        return self._meshed_region.available_named_selections

    @property
    def node_ids(self) -> List[int]:
        """Returns the list of node IDs in the mesh."""
        return self._meshed_region.nodes.scoping.ids

    @property
    def element_ids(self) -> List[int]:
        """Returns the list of element IDs in the mesh."""
        return self._meshed_region.elements.scoping.ids

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.MeshedRegion` object."""
        return self._meshed_region
