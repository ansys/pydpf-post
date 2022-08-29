"""Contains the common ``Result`` class.

This module contains the super class of the stress, strain, structural temperature,
and displacement objects.
"""

from ansys.dpf.core import Operator

from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.result_definition import Definition


class Result:
    """Provides the ``Result`` class.

    This is an abstract class. It is not instantiated directly but instead is
    subclassed by specific result classes like the
    :class:`Displacement <ansys.dpf.post.displacement.Displacement>` class.

    """

    def __init__(self, data_sources, model, **kwargs):
        """Initialize this class."""
        self._data_sources = data_sources
        self._model = model
        self.definition = Definition(**kwargs)

    def __str__(self):
        """Return the string representation of this class."""
        return self.definition.__str__()

    def has_complex_frequencies(self):
        """Check if the result contains complex frequencies.

        Returns
        -------
        bool
            ``True`` if the result contains complex frequencies, ``False`` otherwise.
        """
        tfq = self._model.metadata.time_freq_support
        return tfq.complex_frequencies != None

    def _get_amplitude_evaluation(self, result_data):
        # resultData = self._get_result_data_function_of_operator(
        #     name, self, self._data_sources, **kwargs
        # )
        resultData = result_data
        modulus_op = Operator("modulus")
        modulus_op.inputs.fields_container.connect(
            resultData._evaluator._result_operator.outputs.fields_container
        )
        resultData._evaluator._chained_operators[modulus_op.name] = (
            """This operator computes the amplitude """
            """of the result (when the result has complex values)."""
        )
        # resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        resultData._evaluator._result_operator = modulus_op
        return resultData

    def _get_result_data(
        self, operator_name, data_sources, model, subresult=None, phase=None
    ):
        """Check the keyword arguments that are specified while calling a subresult method.

        The arguments that can be set at this point are:
            - time
            - set
            - phase (if complex result)
        """
        # write correct arguments regarding location
        b_elem_average = False
        location_to_compute = self.definition._location
        # if self.definition._location == locations.elemental_nodal:
        #     location_to_compute = locations.elemental
        # if self.definition._location == locations.elemental:
        #     b_elem_average = True

        return ResultData(
            operator_name=operator_name,
            data_sources=data_sources,
            model=model,
            elem_average=b_elem_average,
            location=location_to_compute,
            element_scoping=self.definition.element_scoping,
            node_scoping=self.definition.node_scoping,
            named_selection=self.definition.named_selection,
            time=self.definition.time,
            grouping=self.definition.grouping,
            subresult=subresult,
            mapdl_grouping=self.definition.mapdl_grouping,
            set=self.definition.set,
            path=self.definition.path,
            time_scoping=self.definition.time_scoping,
            phase=phase,
        )
