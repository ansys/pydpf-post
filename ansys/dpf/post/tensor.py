"""Module containing the Result subclass : Tensor."""

from ansys.dpf.post.result_object import Result

class Tensor(Result):
    """Child class of the Result one.
    Implements a tensor result (stress, strain).
    """
    @property
    def xx(self):
        """Returns XX component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X")

    @property
    def yy(self):
        """Returns YY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y")

    @property
    def zz(self):
        """Returns ZZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z")

    @property
    def xy(self):
        """Returns XY component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XY")

    @property
    def yz(self):
        """Returns YZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="YZ")

    @property
    def xz(self):
        """Returns XZ component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XZ")

    @property
    def principal_1(self):
        """Returns first principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="1")

    @property
    def principal_2(self):
        """Returns second principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="2")

    @property
    def principal_3(self):
        """Returns third principal component of the tensor as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="3")

    @property
    def tensor(self):
        """Returns the tensor values as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model)

    def __str__(self):
        txt = "Tensor object. \n\n"
        txt += super().__str__()
        return txt


class ComplexTensor(Tensor):
    """Child class of the Result one.
    Implements a tensor result (stress,
    strain).
    """
    @property
    def xx_amplitude(self):
        """Returns XX component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X")
        return Result._get_amplitude_evaluation(self, res_data)

    def xx_at_phase(self, phase: float):
        """Returns XX component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="X", phase=phase)

    @property
    def yy_amplitude(self):
        """Returns YY component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y")
        return Result._get_amplitude_evaluation(self, res_data)

    def yy_at_phase(self, phase: float):
        """Returns YY component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Y", phase=phase)

    @property
    def zz_amplitude(self):
        """Returns ZZ component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z")
        return Result._get_amplitude_evaluation(self, res_data)

    def zz_at_phase(self, phase: float):
        """Returns XX component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="Z", phase=phase)

    @property
    def xy_amplitude(self):
        """Returns XY component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XY")
        return Result._get_amplitude_evaluation(self, res_data)

    def xy_at_phase(self, phase: float):
        """Returns XY component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XY", phase=phase)

    @property
    def yz_amplitude(self):
        """Returns YZ component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="YZ")
        return Result._get_amplitude_evaluation(self, res_data)

    def yz_at_phase(self, phase: float):
        """Returns YZ component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="YZ", phase=phase)

    @property
    def xz_amplitude(self):
        """Returns XZ component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XZ")
        return Result._get_amplitude_evaluation(self, res_data)

    def xz_at_phase(self, phase: float):
        """Returns XZ component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="XZ", phase=phase)

    @property
    def principal_1_amplitude(self):
        """Returns first principal component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="1")
        return Result._get_amplitude_evaluation(self, res_data)

    def principal_1_at_phase(self, phase: float):
        """Returns first principal component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="1", phase=phase)

    @property
    def principal_2_amplitude(self):
        """Returns second principal component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="2")
        return Result._get_amplitude_evaluation(self, res_data)

    def principal_2_at_phase(self, phase: float):
        """Returns second principal component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="2", phase=phase)

    @property
    def principal_3_amplitude(self):
        """Returns third principal component of the tensor as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="3")
        return Result._get_amplitude_evaluation(self, res_data)

    def principal_3_at_phase(self, phase: float):
        """Returns third principal component of the tensor at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, subresult="3", phase=phase)

    @property
    def tensor_amplitude(self):
        """Returns the tensor values as a ResultData."""
        res_data = super()._get_result_data(self._operator_name, self._data_sources, self._model)
        return Result._get_amplitude_evaluation(self, res_data)

    def tensor_at_phase(self, phase: float):
        """Returns the tensor values at specific phase as a ResultData."""
        return super()._get_result_data(self._operator_name, self._data_sources, self._model, phase=phase)

    def __str__(self):
        txt = "Complex tensor object. \n\n"
        txt += super().__str__()
        return txt
