"""Module containing the Result subclass : Tensor."""

from ansys.dpf.post.result_object import Result

class Tensor(Result):
    """Child class of the Result one.
    Implements a tensor result (stress, 
    strain).
    """
    @property
    def xx(self):
        """Returns XX component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X")
    
    @property
    def yy(self):
        """Returns YY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y")
    
    @property
    def zz(self):
        """Returns ZZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z")
    
    @property
    def xy(self):
        """Returns XY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XY")
    
    @property
    def yz(self):
        """Returns YZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="YZ")
    
    @property
    def xz(self):
        """Returns XZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XZ")
    
    @property
    def principal_1(self):
        """Returns first principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="1")
    
    @property
    def principal_2(self):
        """Returns second principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="2")
    
    @property
    def principal_3(self):
        """Returns third principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="3")
    
    @property    
    def tensor(self):
        """Returns the tensor values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model)
    
    def __str__(self):
        txt = "Tensor object. \n\n"
        txt += super().__str__()
        return txt