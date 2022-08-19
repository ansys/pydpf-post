"""Module containing the Result subclass : Vector."""

from ansys.dpf.core import Operator

from ansys.dpf.post.result_object import Result


class Vector(Result):
    """Vector (displacement) result."""

    @property
    def x(self):
        """Return X component of the vector as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="X"
        )

    @property
    def y(self):
        """Return Y component of the vector as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Y"
        )

    @property
    def z(self):
        """Return Z component of the vector as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Z"
        )

    @property
    def vector(self):
        """Return the vector values as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )

    @property
    def norm(self):
        """Return the norm of the vector as a ResultData."""
        result_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[
            norm_op.name
        ] = """This operator will compute the norm of the result."""
        return result_data

    def __str__(self):
        """Return the string representation of this class."""
        return f"Vector object.\n\n{super().__str__()}"


class ComplexVector(Vector):
    """Complex Vector."""

    @property
    def x_amplitude(self):
        """Return X component amplitude of the vector as a ResultData."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="X"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def x_at_phase(self, phase: float):
        """Return the X component at specific phase as a ResultData."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="X",
            phase=phase,
        )

    @property
    def y_amplitude(self):
        """Return Y component amplitude of the vector as a ResultData."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Y"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def y_at_phase(self, phase: float):
        """Return the Y component at specific phase as a ResultData."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="Y",
            phase=phase,
        )

    @property
    def z_amplitude(self):
        """Return Z component amplitude of the vector as a ResultData."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Z"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def z_at_phase(self, phase: float):
        """Return the Z component at specific phase as a ResultData."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="Z",
            phase=phase,
        )

    @property
    def vector_amplitude(self):
        """Return the vector amplitude values as a ResultData."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def vector_at_phase(self, phase: float):
        """Return the vector values at specific phase as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )

    @property
    def norm_amplitude(self):
        """Return the amplitude of the norm of the vector as a ResultData."""
        result_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[
            norm_op.name
        ] = """This operator will compute the norm of the result."""
        return Result._get_amplitude_evaluation(self, result_data)

    def norm_at_phase(self, phase: float):
        """Return the norm of the vector at specific phase as a ResultData."""
        result_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[
            norm_op.name
        ] = """This operator will compute the norm of the result."""
        return result_data

    def __str__(self):
        """Return the string representation of this class."""
        return f"Complex vector object.\n\n{super().__str__()}"
