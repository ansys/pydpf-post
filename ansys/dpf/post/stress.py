"""This module contains the stress result classes."""

from ansys.dpf.post.tensor import Tensor, ComplexTensor
from ansys.dpf.post.result_object import Result


class Stress(Tensor):
    """Defines the stress object, that is a tensor object.

    Examples
    --------
    Extract the stress from a solution

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)
    >>> stress = solution.stress(location=post.locations.nodal, time_scoping=[1])
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "S"

    def __str__(self):
        txt = "Stress "
        txt += super().__str__()
        return txt

    @property
    def von_mises(self):
        """Returns the von mises stress as a ResultData."""
        return super()._get_result_data("S_eqv", self._data_sources, self._model)


class ComplexStress(ComplexTensor, Stress):
    """Defines the complex stress object, that is a tensor object."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._operator_name = "S"

    def __str__(self):
        txt = "Complex stress. \n"
        txt += super().__str__()
        return txt

    @property
    def von_mises_amplitude(self):
        """Returns the von mises stress amplitude as a ResultData."""
        res_data = super()._get_result_data("S_eqv", self._data_sources, self._model)
        return Result._get_amplitude_evaluation(self, res_data)

    def von_mises_at_phase(self, phase: float):
        """Returns the von mises stress at specific phase as a ResultData."""
        return super()._get_result_data("S_eqv", self._data_sources, self._model, phase=phase)
