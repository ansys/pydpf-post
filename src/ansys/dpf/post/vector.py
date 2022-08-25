"""Module containing ``Result` subclasses for vectors."""

from ansys.dpf.core import Operator

from ansys.dpf.post.result_object import Result


class Vector(Result):
    """Implements the vector (displacement) result."""

    @property
    def x(self):
        """Result data for the X component of the vector."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="X"
        )

    @property
    def y(self):
        """Result data for the Y component of the vector."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Y"
        )

    @property
    def z(self):
        """Result data for the Z component of the vector."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Z"
        )

    @property
    def vector(self):
        """Result data for the vector values."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )

    @property
    def norm(self):
        """Result data for the norm of the vector."""
        result_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[
            norm_op.name
        ] = """This operator computes the norm of the result."""
        return result_data

    def __str__(self):
        """Return the string representation of this class."""
        return f"Vector object.\n\n{super().__str__()}"


class ComplexVector(Vector):
    """Implements the complex vector."""

    @property
    def x_amplitude(self):
        """Result data for the X component amplitude of the vector."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="X"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def x_at_phase(self, phase: float):
        """Result data for the X component at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="X",
            phase=phase,
        )

    @property
    def y_amplitude(self):
        """Result data for the Y component amplitude of the vector."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Y"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def y_at_phase(self, phase: float):
        """Result data for the Y component at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="Y",
            phase=phase,
        )

    @property
    def z_amplitude(self):
        """Result data for the Z component amplitude of the vector."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Z"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def z_at_phase(self, phase: float):
        """Result data for the Z component at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="Z",
            phase=phase,
        )

    @property
    def vector_amplitude(self):
        """Result data for the vector amplitude values."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def vector_at_phase(self, phase: float):
        """Result data for the vector values at a specific phase."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )

    @property
    def norm_amplitude(self):
        """Result data for the amplitude of the norm of the vector."""
        result_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[
            norm_op.name
        ] = """This operator computes the norm of the result."""
        return Result._get_amplitude_evaluation(self, result_data)

    def norm_at_phase(self, phase: float):
        """Result data for the norm of the vector at a specific phase."""
        result_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )
        out_op = result_data._evaluator._result_operator
        norm_op = Operator("norm_fc")
        norm_op.inputs.fields_container.connect(out_op.outputs.fields_container)
        result_data._evaluator._result_operator = norm_op
        result_data._evaluator._chained_operators[
            norm_op.name
        ] = """This operator computes the norm of the result."""
        return result_data

    def __str__(self):
        """Return the string representation of this class."""
        return f"Complex vector object.\n\n{super().__str__()}"
