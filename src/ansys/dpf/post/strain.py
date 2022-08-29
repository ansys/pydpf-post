"""This module contains classes for strain results."""

from ansys.dpf.post.tensor import ComplexTensor, Tensor


class ElasticStrain(Tensor):
    """Defines the elastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPEL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nElastic strain object."


class ComplexElasticStrain(ComplexTensor):
    """Defines the complex elastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPEL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nComplex elastic strain object."


class PlasticStrain(Tensor):
    """Defines the plastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPPL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nPlastic strain object."


class ComplexPlasticStrain(ComplexTensor):
    """Defines the complex plastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPPL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nComplex plastic strain object."
