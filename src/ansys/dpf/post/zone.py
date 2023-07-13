"""Module containing the ``Zone`` class and sub-classes.

Zone
----

"""
from abc import ABC
from typing import List


class Zone(ABC):
    """Base class for fluid Zone objects."""

    def __init__(
        self, name: str, id: int, scoping
    ):  # pylint: disable=redefined-builtin
        """Initialize this class."""
        self._name = name
        self._id = id
        self._mesh_scoping = scoping

    @property
    def name(self) -> str:
        """Return the name of the Zone."""
        return self._name

    @property
    def id(self) -> int:
        """Return the ID of the Zone."""
        return self._id


class CellZone(Zone):
    """Fluid zone defined on cells."""

    @property
    def cells(self) -> List[int]:
        """Return the list of IDs of the cells in the Zone."""
        return self._mesh_scoping.ids


class FaceZone(Zone):
    """Fluid zone defined on faces."""

    @property
    def faces(self) -> List[int]:
        """Return the list of IDs of the faces in the Zone."""
        return self._mesh_scoping.ids


class Zones:
    """List of fluid zones."""

    def __init__(self):
        """Initialize this class."""
        self._zones = []  # list of CellZone and FaceZone
        self._face_zones_ind = 0

    @property
    def cell_zones(self) -> List[CellZone]:
        """Return the list of CellZone instances only."""
        return self._zones[: self._face_zones_ind]

    def face_zones(self) -> List[FaceZone]:
        """Return the list of FaceZone instances only."""
        return self._zones[self._face_zones_ind :]

    def __getitem__(self, key: int) -> Zone:
        """Return the Zone at the given position in the list."""
        return self._zones[key]
