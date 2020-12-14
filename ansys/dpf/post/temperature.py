"""This module contains the temperature result class ."""

from ansys.dpf.post.scalar import Scalar, ComplexScalar


class Temperature(Scalar):
    """Defines the strctural temperature object, that is a scalar object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a temperature object."
        return txt
    
    
class ComplexTemperature(ComplexScalar):
    """Defines the complex strctural temperature object, that is a scalar object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a temperature object."
        return txt