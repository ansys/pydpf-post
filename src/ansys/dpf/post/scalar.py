"""Module containing the Result subclass : Scalar."""

from ansys.dpf.post.result_object import Result


class Scalar(Result):
    """Child class of the Result.

    Implements a scalar result (temperature).
    """

    @property
    def scalar(self):
        """Return the scalar values as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )

    def __str__(self):
        """Return the string representation."""
        return f"Scalar object.\n\n{super().__str__()}"


class ComplexScalar(Scalar):
    """Child class of the Result.

    Implements a complex scalar result (temperature).
    """

    @property
    def scalar_amplitude(self):
        """Return the scalar amplitude values as a ResultData."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def scalar_at_phase(self, phase: float):
        """Return the scalar values at specified phase as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )

    def __str__(self):
        """Return the string representation."""
        return f"Complex scalar object.\n\n{super().__str__()}"
