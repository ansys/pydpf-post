"""This module contains the super class of the
stress/strain/structural_temperature/displacement objects."""

from ansys.dpf.core.common import locations
from ansys.dpf.core import Operator
from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.result_definition import Definition

class Result:
    def __init__(self, data_sources, model, **kwargs):
        self._data_sources = data_sources
        self._model = model
        self.definition = Definition(**kwargs)
        self._op_average = None

    def __str__(self):
        txt = self.definition.__str__()
        return txt

    def has_complex_frequencies(self):
        """Returns True if the result support contains complex frequencies."""
        tfq = self._model.metadata.time_freq_support
        if (tfq.complex_frequencies != None):
            return True
        else:
            return False

    def _get_amplitude_evaluation(self, result_data):
        # resultData = self._get_result_data_function_of_operator(name, self, self._data_sources, **kwargs)
        resultData = result_data
        modulus_op = Operator("modulus")
        modulus_op.inputs.fields_container.connect(resultData._evaluator._result_operator.outputs.fields_container)
        resultData._evaluator._chained_operators[modulus_op.name] = """This operator will compute the amplitude of the result (when result has complex values)."""
        # resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        resultData._evaluator._result_operator = modulus_op
        return resultData

    def _get_result_data(self, operator_name, data_sources, model, subresult=None, phase=None):
        """This method checks the keyword arguments that are
        specified while calling a subresult method.

        The arguments that can be set at this point are:
            - time
            - set
            - phase (if complex result)
        """
        #write correct arguments regarding location
        b_elem_average = False
        location_to_compute = self.definition._location
        if self.definition._location == locations.elemental_nodal:
            location_to_compute = locations.elemental
        if self.definition._location == locations.elemental:
            b_elem_average = True

        return ResultData(operator_name=operator_name, data_sources=data_sources,
                          model=model, elem_average=b_elem_average, op_average=self._op_average, location=location_to_compute,
                          element_scoping=self.definition.element_scoping, node_scoping=self.definition.node_scoping,
                          named_selection=self.definition.named_selection,
                          time=self.definition.time, grouping=self.definition.grouping,
                          subresult=subresult, mapdl_grouping=self.definition.mapdl_grouping, set=self.definition.set,
                          time_scoping=self.definition.time_scoping, phase=phase)
