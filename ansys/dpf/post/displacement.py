"""This module contains the displacement result class ."""

from ansys.dpf.post.vector import Vector


class Displacement(Vector):
    """Defines the displacement object, that is a vector object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "U"
        
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a displacement result."
        return txt