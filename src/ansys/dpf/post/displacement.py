"""This module contains the displacement result class."""

from ansys.dpf.core import locations

from ansys.dpf.post.errors import NodalLocationError
from ansys.dpf.post.vector import ComplexVector, Vector


class Displacement(Vector):
    """Defines the displacement object, which is a vector object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "U"

        if self.definition.location != locations.nodal:
            raise NodalLocationError

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nDisplacement object."


class ComplexDisplacement(ComplexVector):
    """Defines the complex displacement object, which is a vector object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "U"

        if self.definition.location != locations.nodal:
            raise NodalLocationError
