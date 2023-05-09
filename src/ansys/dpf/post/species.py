"""Module containing the ``Species`` class and sub-classes.

Species
-------
In fluid simulations, the species are the chemical compounds present.
"""


class Species:
    """In fluid simulation, defines a chemical species."""

    def __init__(self, name: str, id: int):
        """Initialize this class."""
        self._name = name
        self._id = id

    @property
    def name(self) -> str:
        """Return the name of the species."""
        return self._name

    @property
    def id(self) -> int:
        """Return the ID of the species."""
        return self._id

    def __repr__(self) -> str:
        """String representation of the instance."""
        return f'Species<name: "{self._name}", id={self._id}>'


class SpeciesList:
    """In fluid simulation, defines a list of chemical species."""

    def __init__(self):
        """Initialize this class."""
        self._species = []

    def __repr__(self) -> str:
        """String representation of the instance."""
        return "[".join([repr(species) for species in self]) + "]"

    def __getitem__(self, item: int) -> Species:
        """Return the Species at the given position in the list."""
        return self._species[item]
