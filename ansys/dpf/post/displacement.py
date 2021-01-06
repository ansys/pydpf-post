"""This module contains the displacement result class."""

from ansys.dpf.post.vector import Vector, ComplexVector
from ansys.dpf.core import locations
from ansys.dpf.post.errors import NodalLocationError

class Displacement(Vector):
    """Defines the displacement object, that is a vector object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "U"
        
        if self.definition.location != locations.nodal:
            raise NodalLocationError
        
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a displacement object."
        return txt
    
    
class ComplexDisplacement(ComplexVector):
    """Defines the complex displacement object, that is a vector object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "U"
        
        if self.definition.location != locations.nodal:
            raise NodalLocationError
