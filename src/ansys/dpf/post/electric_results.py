"""This module contains the electric results class ."""

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.post.scalar import Scalar
from ansys.dpf.post.vector import Vector


class ElectricField(Vector):
    """Defines the temperature object for a thermal/electric analysis."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EF"

        # disable element scoping
        if _AvailableKeywords.element_scoping in kwargs:
            raise Exception(
                "Element scoping is not available with thermal/electric results."
            )
        self.definition._Definition__element_scoping_locked = True

    def __str__(self):
        """Return the string representation."""
        return f"{super().__str__()}\nElectric field object."


class ElectricPotential(Scalar):
    """Defines the temperature object for a thermal/electric analysis."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "VOLT"

        # disable element scoping
        if _AvailableKeywords.element_scoping in kwargs:
            raise Exception(
                "Element scoping is not available with thermal/electric results."
            )
        self.definition._Definition__element_scoping_locked = True

    def __str__(self):
        """Return the string representation."""
        return f"{super().__str__()}\nElectric potential object."
