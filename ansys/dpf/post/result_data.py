"""Module containing the _ResultData class ("result object") that
user will be able to use to compute through the DPF Post API.

This is a fields container wrapper."""

from textwrap import wrap
import os
import sys

import pyvista as pv
import matplotlib.pyplot as pyplot

from ansys.dpf.core import Operator, FieldsContainer
from ansys.dpf.core.common import types
from ansys.dpf.core.plotter import Plotter as DpfPlotter
from ansys.dpf.post.result_evaluation import ResultEvaluator


class ResultData:
    """Instance of the result of a dpf result object.
    Created thanks to the dpf result object.

    Parameters
    ----------
    Following list of keywords:
        - location
        - node_scoping
        - element_scoping
        - time
        - grouping

    The entire list of parameters can be found using
     ``post.print_available_keywords()``.

    Examples
    --------
    >>> from ansys.dpf import post
    >>> solution = post.solution("file.rst")
    >>> disp = solution.nodal_displacement()
    >>> disp_on_nodes = solution.nodal_displacement(node_scoping = [1, 23])
    >>> disp_on_named_selection = solution.nodal_displacement(named_selection="SELECTION")
    """
    
    def __init__(self, operator_name: str, data_sources, model, 
                 elem_average: bool, op_average = None, 
                 location: str = None, element_scoping = None, 
                 node_scoping = None, named_selection = None, 
                 time = None, 
                 grouping = None, phase = None, subresult = None, 
                 mapdl_grouping = None, set = None, time_scoping = None):
        
        self._evaluator = ResultEvaluator(operator_name, data_sources, model, elem_average, 
                 op_average, location, element_scoping, node_scoping, named_selection, time, 
                 grouping, phase, subresult, mapdl_grouping, set, time_scoping)
        self.result_fields_container = None
             
            
    def __str__(self):
        self._evaluate_result()
        name = self.result_fields_container[0].name.split("_")
        txt = "%s result.\n" % name[0].capitalize()
        if (self._evaluator.subresult is not None):
            txt += "%s component. \n" % self._evaluator.subresult
        txt += "\n"
        desc = "This result has been computed using dpf.core.Operator objects, which "
        desc += "have been chained together according to the following list:"
        txt += '\n'.join(wrap(desc)) + '\n'
        for key, val in self._evaluator._chained_operators.items():
            txt += "- %s: " % key
            txt += val
            txt += "\n"
        # txt += "\n\n"
        # txt += self._evaluator._model.__str__()
        return txt

    def _evaluate_result(self):
        """First evaluation of the result."""
        if self.result_fields_container is None:
            self.result_fields_container = self._evaluator.evaluate_result()

    def _evaluate_result_forced(self):
        """Re-evaluation of the result."""
        self.result_fields_container = self._evaluator.evaluate_result()

    def get_all_label_spaces(self):
        """Returns all the label spaces contained in a result 
        as a string.
        Labels can be used to select fields to plot.
        
        Returns
        -------
        list 
            List of dictionary (list of label space)
        """
        self._evaluate_result()
        i = 0
        list_labels = []
        while i < len(self.result_fields_container):
            list_labels.append(self.result_fields_container.get_label_space(i))
            i += 1
        return list_labels
    
    @property
    def num_fields(self):
        """Returns the number of fields contained in the result."""
        self._evaluate_result()
        return len(self.result_fields_container)
    
    def get_data_at_field(self, field_index: int = 0):
        """Returns the data at the field with the specified index."""
        self._evaluate_result()
        return self.result_fields_container[field_index].data
    
    def __getitem__(self, field_index: int = 0):
        """Override of the result item getter. Implements the fields_container_result 
        item getter. Return the field at the field_index position in the result fields 
        container.
        """
        self._evaluate_result()
        return self.result_fields_container[field_index]
    
    def get_scoping_at_field(self, field_index: int = 0):
        """Returns the scoping of the result."""
        self._evaluate_result()
        field = self.result_fields_container[field_index]
        return field.scoping.ids
        
    def _min_max(self, pin):
        """"Chains the operators to compute the min/max values."""
        self._evaluate_result()
        max_operator = Operator("min_max_fc")
        max_operator.inputs.connect(self.result_fields_container)
        result = max_operator.get_output(pin, types.field)
        return result
    
    @property    
    def max(self):
        """Returns the maximum values field. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 1).
        """
        return self._min_max(1)
    
    @property
    def max_data(self):
        """Returns the maximum values field data. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 1).
        """
        return self._min_max(1).data
    
    def get_max_data_at_field(self, field_index: int = 0):
        """Returns the maximum values field data at field_index. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 1).
        """
        return self._min_max(1).data[field_index]
    
    @property
    def min(self):
        """Returns the minimum values field. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 0).
        """
        return self._min_max(0)
    
    @property
    def min_data(self):
        """Returns the minimum values field data. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 0).
        """
        return self._min_max(0).data
    
    def get_min_data_at_field(self, field_index: int = 0):
        """Returns the minimum values field data at field_index. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 0).
        """
        return self._min_max(0).data[field_index]
    
    
    def _plot_contour_with_vtk_file(self):
        """Plot the contour result on its mesh support. The obtained figure depends on the 
        support (can be a meshed_region or a time_freq_support).
        If transient analysis, plot the last result.
        
        This method is private, publishes a vtk file and print (using pyvista) from this file."""
        self._evaluate_result()
        pl = DpfPlotter(self._evaluator._model.metadata.meshed_region)
        pl._plot_contour_using_vtk_file(self.result_fields_container)
        
    
    def plot_contour(self, display_option: str = "time", option_id=1,
                     off_screen=None, notebook=None, **kwargs):
        """Plot the contour result on its mesh support.

        The obtained figure depends on the support (can be a
        meshed_region or a time_freq_support).  If transient analysis,
        plot the last result if no time_scoping has been specified.
        The self.get_all_label_spaces() will return a string
        containing all the label spaces.

        Parameters
        ----------
        display_option : str, optional
            The name of the label you want to display. Default is ``"time"``.
        option_id: int, optional
             The list of label ids you want to display. Default is ``[1]``.
        off_screen : bool, optional
            Renders off screen when ``True``.  Useful for automated screenshots.
        notebook : bool, optional
            Option to force plotting within a jupyter notebook or
            external to the notebook as an interactive figure.
        **kwargs : optional
            Additional keyword arguments for the plotter.  See
            ``help(pyvista.plot)`` for additional keyword arguments.

        Examples
        --------
        Plot a result at the time_step number 1

        >>> from ansys.dpf import post
        >>> solution = post.load_solution('file.rst')
        >>> stress = solution.stress(location=post.locations.nodal)
        >>> sx = stress.xx
        >>> sx.plot_contour("time", [1])

        The labels can be obtained using:

        >>> sx.get_all_label_spaces()

        [{'elshape': 1, 'time': 1}, {'elshape': 0, 'time': 1}]
        """
        self._evaluate_result()
        pl = DpfPlotter(self._evaluator._model.metadata.meshed_region)
        if len(self.result_fields_container) == 1:
            pl.plot_contour(self.result_fields_container, off_screen=off_screen,
                            notebook=notebook, **kwargs)
        else:
            # sorts and creates a new fields_container with only the desired labels
            ids = option_id
            if isinstance(option_id, int):
                ids = [option_id]
            new_fields_container = FieldsContainer()
            label_spaces = []
            for opt_id in ids:
                i = 0
                while i < len(self.result_fields_container):
                    label_space = {}
                    label_space[display_option] = opt_id
                    fc_label = self.result_fields_container.get_label_space(i)
                    for lab in self.result_fields_container.labels:
                        if lab != display_option:
                            if display_option in fc_label:
                                if fc_label[display_option] == opt_id:
                                    label_space[lab] = fc_label[lab]
                    label_spaces.append(label_space)
                    i += 1
            lab_names = []
            for lab_sp in label_spaces:
                for key, val in lab_sp.items():
                    if key not in lab_names:
                        lab_names.append(key)
            for name in lab_names:
                new_fields_container.add_label(name)
            for label in label_spaces:
                field = self.result_fields_container._get_entries(label)
                new_fields_container.add_field(label, field)
            pl.plot_contour(new_fields_container, off_screen=off_screen,
                            notebook=notebook, **kwargs)

    def _plot_chart(self):
        """Plot the minimum/maximum result values over time.

        This method works if the time_freq_support contains
        several time_steps (for example, a transient analysis).

        A time_scoping keyword must be used to select all the
        time_steps of the result.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution('file.rst')
        >>> tscope = list(range(1, len(solution.time_freq_support.frequencies) + 1))
        >>> stress = solution.stress(mapdl_grouping=181, location='Nodal',
                                     time_scoping=tscope
        >>> s = stress.tensor
        >>> s.plot_chart()
        """
        self._evaluate_result()
        # tfq = self._evaluator._model.metadata.time_freq_support
        # timeids = list(range(1,tfq.n_sets+1))
        # res_op = self._evaluator._result_operator
        # res_op.inputs.time_scoping.connect(timeids)
        # new_fields_container = res_op.get_output(0, types.fields_container)
        pl = DpfPlotter(self._evaluator._model.metadata.meshed_region)
        pl.plot_chart(self.result_fields_container)
