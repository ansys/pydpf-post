"""Module containing the Result subclass : Scalar."""

from ansys.dpf.post.result_object import Result

class Scalar(Result):
    """Child class of the Result one.
    Implements a scalar result (temperature).
    """
    
    @property
    def scalar(self):
        """Returns the scalar values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model)
    
    def __str__(self):
        txt = "Scalar object. \n\n"
        txt += super().__str__()
        return txt
    
    
class ComplexScalar(Scalar):
    """Child class of the Result one.
    Implements a complex scalar result (temperature).
    """
    
    @property
    def scalar_amplitude(self):
        """Returns the scalar amplitude values as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model)
        return Result._get_amplitude_evaluation(self, res_data)
    
    def scalar_at_phase(self, phase: float):
        """Returns the scalar values at specified phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, phase=phase)
    
    def __str__(self):
        txt = "Complex scalar object. \n\n"
        txt += super().__str__()
        return txt