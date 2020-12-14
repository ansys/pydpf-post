"""This module contains the strain result class ."""

from ansys.dpf.post.tensor import Tensor, ComplexTensor


class ElasticStrain(Tensor):
    """Defines the elastic strain object, that is a tensor object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "EPEL"
        
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is an elastic strain object."
        return txt
    

class ComplexElasticStrain(ComplexTensor):
    """Defines the complex elastic strain object, that is a tensor object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "EPEL"
        
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is an elastic strain object."
        return txt
        
        
class PlasticStrain(Tensor):
    """Defines the plastic strain object, that is a tensor object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "EPPL"
        
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a plastic strain object."
        return txt
    

class ComplexPlasticStrain(ComplexTensor):
    """Defines the complex plastic strain object, that is a tensor object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "EPPL"
        
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This is a plastic strain object."
        return txt