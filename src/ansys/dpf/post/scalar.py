"""Module containing the Result subclass : Scalar."""

from ansys.dpf.post.result_object import Result


class Scalar(Result):
    """Provides a child ``Result`` class that implements scalar results (temperatures)."""

    @property
    def scalar(self):
        """Result data for the scalar values."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )

    def __str__(self):
        """Return the string representation."""
        return f"Scalar object.\n\n{super().__str__()}"


class ComplexScalar(Scalar):
    """Provides a child ``Result`` class that implements complex scalar results (temperatures)."""

    @property
    def scalar_amplitude(self):
        """Get result data for the scalar amplitude values."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def scalar_at_phase(self, phase: float):
        """Get result data for the scalar values at a specified phase."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )

    def __str__(self):
        """Return the string representation."""
        return f"Complex scalar object.\n\n{super().__str__()}"
