"""Module containing the Result subclass : Vector."""

from ansys.dpf.core import Operator
from ansys.dpf.post.result_object import Result

class Vector(Result):
    """Child class of the Result one.
    Implements a vector result (displacement).
    """
    @property
    def x(self):
        """Returns X component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X")
    
    @property
    def y(self):
        """Returns Y component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y")
    
    @property
    def z(self):
        """Returns Z component of the vector as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z")
    
    @property
    def vector(self):
        """Returns the vector values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model)
    
    @property
    def norm(self):
        """Returns the norm of the vector as a ResultData."""
        result_data = super()._get_result_data(self._operator_name, self._data_sources, self._model)    
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