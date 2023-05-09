"""Module containing the ``Phase`` class and sub-classes.

Phase
-----
In fluid simulations, a phase defines the physical state of a given species.
"""


class Phase:
    """Physical state of a given species in a fluid simulation."""

    def __init__(self, name: str, id: int):
        """Initialize this class."""
        self._name = name
        self._id = id

    @property
    def name(self) -> str:
        """Return the name of the phase."""
        return self._name

    @property
    def id(self) -> int:
        """Return the ID of the phase."""
        return self._id

    def __repr__(self) -> str:
        """String representation of the instance."""
        return f'Phase<name: "{self._name}", id={self._id}>'


class Phases:
    """List of physical states in a fluid simulation."""

    def __init__(self):
        """Initialize this class."""
        self._phases = []

    def __repr__(self) -> str:
        """String representation of the instance."""
        return "[".join([repr(phase) for phase in self]) + "]"

    def __getitem__(self, item: int) -> Phase:
        """Return the Phase at the given position in the list."""
        return self._phases[item]
