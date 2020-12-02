"""Module containing the Result subclass : Vector."""

from ansys.dpf.post.result_object import Result

class Vector(Result):
    """Child class of the Result one.
    Implements a vector result (displacement).
    """
    def X(self, **kwargs):
        """Returns X component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X", **kwargs)
    
    def Y(self, **kwargs):
        """Returns Y component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y", **kwargs)
    
    def Z(self, **kwargs):
        """Returns Z component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z", **kwargs)
    
    def vector(self, **kwargs):
        """Returns the vector values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, **kwargs)
    
    def __str__(self):
        txt = "Vector object. \n\n"
        txt += super().__str__()
        return txt