"""This module contains classes for temperature results ."""

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.post.scalar import ComplexScalar, Scalar
from ansys.dpf.post.vector import Vector


class StructuralTemperature(Scalar):
    """Defines the structural temperature object, which is a scalar object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nStructural temperature object."


class ComplexStructuralTemperature(ComplexScalar):
    """Defines the complex structural temperature object, which is a scalar object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nComplex temperature object."


class Temperature(Scalar):
    """Defines the temperature object, which is a scalar object, for thermal analysis."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "TEMP"

        # disable element scoping
        if _AvailableKeywords.element_scoping in kwargs:
            raise Exception(
                "Element scoping is not available with thermal/electric results."
            )
        self.definition._Definition__element_scoping_locked = True

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nTemperature object."


class HeatFlux(Vector):
    """Defines the heat flux object, which is a scalar object, for thermal analysis."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "TF"

        # disable element scoping
        if _AvailableKeywords.element_scoping in kwargs:
            raise Exception(
                "Element scoping is not available with thermal/electric results."
            )
        self.definition._Definition__element_scoping_locked = True

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nHeat flux object."
