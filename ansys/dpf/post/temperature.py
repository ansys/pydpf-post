"""This module contains the temperature result class ."""

from ansys.dpf.post.scalar import Scalar, ComplexScalar


class StructuralTemperature(Scalar):
    """Defines the strctural temperature object, that is a scalar object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a temperature object."
        return txt
    
    
class ComplexStructuralTemperature(ComplexScalar):
    """Defines the complex strctural temperature object, that is a scalar object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a temperature object."
        return txt
    
    
class Temperature(Scalar):
    """Defines the temperature object for thermal analysis, that is a scalar object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "TEMP"

    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a temperature object."
        return txt