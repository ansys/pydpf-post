"""This module contains the stress result class ."""

from ansys.dpf.post.tensor import Tensor


class Stress(Tensor):
    """Defines the stress object, that is a tensor object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "S"
        
    def __str__(self):
        txt = "Stress. \n"
        txt += super().__str__()
        return txt
    
    def von_mises(**kwargs):
        pass