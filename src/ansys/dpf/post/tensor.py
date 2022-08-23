"""Module containing the Tensor Result subclass."""

from ansys.dpf.post.result_object import Result


class Tensor(Result):
    """Provides a child ``Result`` class that implements a tensor result (stress, strain)."""

    @property
    def xx(self):
        """Return XX component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="X"
        )

    @property
    def yy(self):
        """Return YY component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Y"
        )

    @property
    def zz(self):
        """Return ZZ component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Z"
        )

    @property
    def xy(self):
        """Return XY component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="XY"
        )

    @property
    def yz(self):
        """Return YZ component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="YZ"
        )

    @property
    def xz(self):
        """Return XZ component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="XZ"
        )

    @property
    def principal_1(self):
        """Return first principal component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="1"
        )

    @property
    def principal_2(self):
        """Return second principal component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="2"
        )

    @property
    def principal_3(self):
        """Return third principal component of the tensor as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="3"
        )

    @property
    def tensor(self):
        """Return the tensor values as a ResultData."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )

    def __str__(self):
        """Return the string representation of this class."""
        return f"Tensor object.\n\n{super().__str__()}"


class ComplexTensor(Tensor):
    """Child ``Result`` class that implements a complex tensor result (stress, strain)."""

    @property
    def xx_amplitude(self):
        """Get the result data for the XX component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="X"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def xx_at_phase(self, phase: float):
        """Get the result data for the XX component of the tensor at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="X",
            phase=phase,
        )

    @property
    def yy_amplitude(self):
        """Get the result data for the YY component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Y"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def yy_at_phase(self, phase: float):
        """Get the result data for the YY component of the tensor at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="Y",
            phase=phase,
        )

    @property
    def zz_amplitude(self):
        """Get the result data for the ZZ component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="Z"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def zz_at_phase(self, phase: float):
        """Get the result data for the XX component of the tensor at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="Z",
            phase=phase,
        )

    @property
    def xy_amplitude(self):
        """Get the result data for the XY component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="XY"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def xy_at_phase(self, phase: float):
        """Get the result data for the XY component of the tensor at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="XY",
            phase=phase,
        )

    @property
    def yz_amplitude(self):
        """Get the result data for the YZ component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="YZ"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def yz_at_phase(self, phase: float):
        """Get the result data for the YZ component of the tensor at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="YZ",
            phase=phase,
        )

    @property
    def xz_amplitude(self):
        """Get the result data for the  XZ component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="XZ"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def xz_at_phase(self, phase: float):
        """Get the result data for the XZ component of the tensor at a specific phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="XZ",
            phase=phase,
        )

    @property
    def principal_1_amplitude(self):
        """Get the result data for the first principal component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="1"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def principal_1_at_phase(self, phase: float):
        """Get the result data for the first principal component of the tensor at a phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="1",
            phase=phase,
        )

    @property
    def principal_2_amplitude(self):
        """Get the result data for the second principal component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="2",
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def principal_2_at_phase(self, phase: float):
        """Get the result data for the second principal component of the tensor at a phase."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="2",
            phase=phase,
        )

    @property
    def principal_3_amplitude(self):
        """Get the result data for the third principal component of the tensor."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model, subresult="3"
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def principal_3_at_phase(self, phase: float):
        """Return third principal component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(
            self._operator_name,
            self._data_sources,
            self._model,
            subresult="3",
            phase=phase,
        )

    @property
    def tensor_amplitude(self):
        """Get the result data for the the tensor values."""
        res_data = super()._get_result_data(
            self._operator_name, self._data_sources, self._model
        )
        return Result._get_amplitude_evaluation(self, res_data)

    def tensor_at_phase(self, phase: float):
        """Get the result data for the the tensor values at a specific phase."""
        return super()._get_result_data(
            self._operator_name, self._data_sources, self._model, phase=phase
        )

    def __str__(self):
        """Return the string representation of this class."""
        return f"Complex tensor object.\n\n{super().__str__()}"
