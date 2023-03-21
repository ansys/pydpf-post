"""Module containing the ``Meshes`` class.

Meshes
------

"""
from __future__ import annotations

from ansys.dpf.core import MeshesContainer

from ansys.dpf.post.mesh import Mesh


class Meshes:
    """Container to hold and interact with a split mesh."""

    def __init__(self, meshes_container: MeshesContainer):
        """Initialize this class."""
        self._core_object = meshes_container

    def __getitem__(self, item) -> Mesh:
        """Select a specific mesh based on its position in the container."""
        return Mesh(meshed_region=self._core_object[item])

    def __str__(self):
        """String representation of this class."""
        return str(self._core_object)
