"""Module containing the Result subclass : Scalar."""

from ansys.dpf.post.result_object import Result

class Scalar(Result):
    """Child class of the Result one.
    Implements a scalar result (temperature).
    """
    def scalar(self, **kwargs):
        """Returns the scalar values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, **kwargs)
    
    def __str__(self):
        txt = "Scalar object. \n\n"
        txt += super().__str__()
        return txt