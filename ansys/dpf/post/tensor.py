"""Module containing the Result subclass : Tensor."""

from ansys.dpf.post.result_object import Result

class Tensor(Result):
    """Child class of the Result one.
    Implements a tensor result (stress, 
    strain).
    """
    def xx(self, **kwargs):
        """Returns XX component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X", **kwargs)
    
    def yy(self, **kwargs):
        """Returns YY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y", **kwargs)
    
    def zz(self, **kwargs):
        """Returns ZZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z", **kwargs)
    
    def xy(self, **kwargs):
        """Returns XY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XY", **kwargs)
    
    def yx(self, **kwargs):
        """Returns YZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="YZ", **kwargs)
    
    def xz(self, **kwargs):
        """Returns XZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XZ", **kwargs)
    
    def principal_1(self, **kwargs):
        """Returns first principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="1", **kwargs)
    
    def principal_2(self, **kwargs):
        """Returns second principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="2", **kwargs)
    
    def principal_3(self, **kwargs):
        """Returns third principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="3", **kwargs)
        
    def tensor(self, **kwargs):
        """Returns the tensor values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, **kwargs)
    
    def __str__(self):
        txt = "Tensor object. \n\n"
        txt += super().__str__()
        return txt