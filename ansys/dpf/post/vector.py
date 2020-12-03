"""Module containing the Result subclass : Vector."""

from ansys.dpf.core import Operator
from ansys.dpf.post.result_object import Result

class Vector(Result):
    """Child class of the Result one.
    Implements a vector result (displacement).
    """
    def x(self, **kwargs):
        """Returns X component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X", **kwargs)
    
    def y(self, **kwargs):
        """Returns Y component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y", **kwargs)
    
    def z(self, **kwargs):
        """Returns Z component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z", **kwargs)
    
    def vector(self, **kwargs):
        """Returns the vector values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, **kwargs)
    
    def norm(self, **kwargs):
        """Returns the norm of the vector as a ResultData."""
        result_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, **kwargs)    
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[norm_op.name] = """This operator will compute the norm of the result."""
        return result_data
    
    def __str__(self):
        txt = "Vector object. \n\n"
        txt += super().__str__()
        return txt