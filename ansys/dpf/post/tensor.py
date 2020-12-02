"""Module containing the Result subclass : Tensor."""

from ansys.dpf.post.result_object import Result

class Tensor(Result):
    """Child class of the Result one.
    Implements a tensor result (stress, 
    strain).
    """
    def XX(self, **kwargs):
        """Returns XX component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X", **kwargs)
    
    def YY(self, **kwargs):
        """Returns YY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y", **kwargs)
    
    def ZZ(self, **kwargs):
        """Returns ZZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z", **kwargs)
    
    def XY(self, **kwargs):
        """Returns XY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XY", **kwargs)
    
    def YZ(self, **kwargs):
        """Returns YZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="YZ", **kwargs)
    
    def XZ(self, **kwargs):
        """Returns XZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XZ", **kwargs)
    
    def tensor(self, **kwargs):
        """Returns the tensor values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, **kwargs)
    
    def __str__(self):
        txt = "Tensor object. \n\n"
        txt += super().__str__()
        return txt