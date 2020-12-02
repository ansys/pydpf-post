"""Module containing the _ResultData class ("result object") that
user will be able to use to compute through the DPF Post API.

This is a fields container wrapper."""

import pyvista as pv
import matplotlib.pyplot as pyplot
import os
import sys
from ansys.dpf.core import Operator
from ansys.dpf.core.common import types
from ansys.dpf.core.rescoper import Rescoper as _Rescoper
from ansys.dpf.post.result_evaluation import ResultEvaluator


class ResultData:
    """Instance of the result of a dpf result object.
    Created thanks to the dpf result object, itself instantiated from post.
    
    Parameters
    -----
    Following list of keywords:
        - location
        - node_scoping
        - element_scoping
        - named_selection
        - phase
        - time_step
        - grouping
        (...)
        
    The whole list of parameters can be found using print(post.available_keywords()).
    
    Example
    -----
    from ansys.dpf import post
    result = post.result("file.rst")
    disp = result.nodal_displacement() 
    disp_on_nodes = result.nodal_displacement(node_scoping = [1, 23]) 
    disp_on_named_selection = result.nodal_displacement(named_selection = "SELECTION")
    """
    
    def __init__(self, operator_name: str, data_sources, model, 
                 elem_average: bool, 
                 location: str = None, element_scoping = None, 
                 node_scoping = None, named_selection = None, 
                 time = None, 
                 grouping = None, phase = None, subresult = None, 
                 mapdl_grouping = None, set = None, time_scoping = None):
        
        self._evaluator = ResultEvaluator(operator_name, data_sources, model, elem_average, 
                 location, element_scoping, node_scoping, named_selection, time, 
                 grouping, phase, subresult, mapdl_grouping, set, time_scoping)
        self.result_fields_container = None
             
            
    def __str__(self):
        self._evaluate_result()
        name = self.result_fields_container[0].name.split("_")
        txt = "%s result.\n" % name[0].capitalize()
        if (self._evaluator.subresult is not None):
            txt += "%s component. \n" % self._evaluator.subresult
        txt += "\n"
        txt += "The result is computed thanks to dpf.core.Operator objects, "
        txt += "that are chained together regarding the following list: \n"
        for key, val in self._evaluator._chained_operators.items():
            txt += "- %s: " % key
            txt += val
            txt += "\n"
        # txt += "\n\n"
        # txt += self._evaluator._model.__str__()
        return txt
    
    
    def _evaluate_result(self):
        """First evaluation of the result."""
        if (self.result_fields_container == None):
            self.result_fields_container = self._evaluator.evaluate_result()
            
    def _evaluate_result_forced(self):
        """Re-evaluation of the result."""
        self.result_fields_container = self._evaluator.evaluate_result()
        
        
    def get_all_label_spaces(self):
        """Returns all the label spaces contained in a result 
        as a string.
        Labels can be used to select fields to plot.
        
        Returns
        -----
        str
        """
        self._evaluate_result()
        i = 0
        txt = ""
        while i < self.result_fields_container.__len__():
            txt += self.result_fields_container.get_label_space(i).__str__()
            txt += "\n"
            i += 1
        return txt
    
    def num_fields(self):
        """Returns the number of fields contained in the result."""
        self._evaluate_result()
        return self.result_fields_container.__len__()
    
    def data_at_field(self, field_index: int = 0):
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
    
    def scoping_at_field(self, field_index: int = 0):
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
        
    def max(self):
        """Returns the maximum values field. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 1).
        """
        return self._min_max(1)
    
    def max_data(self):
        """Returns the maximum values field data. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 1).
        """
        return self._min_max(1).data
    
    def max_data_at_field(self, field_index: int = 0):
        """Returns the maximum values field data at field_index. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 1).
        """
        return self._min_max(1).data[field_index]
    
    def min(self):
        """Returns the minimum values field. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 0).
        """
        return self._min_max(0)
    
    def min_data(self):
        """Returns the minimum values field data. 
        Chains the result operator to the "min_max_fc" operator, 
        returns its result (output from pin 0).
        """
        return self._min_max(0).data
    
    def min_data_at_field(self, field_index: int = 0):
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
        plotter = pv.Plotter()
        mesh_provider = Operator("MeshProvider")
        mesh_provider.inputs.data_sources.connect(self._evaluator._model.metadata.data_sources)
        vtk_export = Operator("vtk_export")
        path = os.getcwd()
        file_name = "dpf_temporary_hokflb2j9sjd0a3.vtk"
        path += "/" + file_name
        vtk_export.inputs.mesh.connect(mesh_provider.outputs.mesh)
        vtk_export.inputs.fields1.connect(self.result_fields_container)
        vtk_export.inputs.file_path.connect(path)
        vtk_export.run()
        grid = pv.read(path)
        if os.path.exists(path):
            os.remove(path)
        names = grid.array_names
        field_name = self.result_fields_container[0].name
        for n in names: #get new name (for example if time_steps)
            if field_name in n:
                field_name = n #default: will plot the last time_step 
        val = grid.get_array(field_name)
        plotter.add_mesh(grid, scalars=val, stitle = field_name, show_edges=True)
        plotter.add_axes()
        plotter.show()
        
    
    def plot_contour(self, display_option: str = "time", option_id: [int] = [1]):
        """Plot the contour result on its mesh support. The obtained figure depends on the 
        support (can be a meshed_region or a time_freq_support).
        If transient analysis, plot the last result if no time_scoping has been specified.
        
        Parameters
        -----
        display_option: str (the name of the label you want to display). Default is "time".
        option_id: int (the list of label ids you want to display). Default is [1].
        
        Help
        -----
        The self.get_all_label_spaces() will return a string containing all the label spaces.
        
        Example
        -----
        The following labels are obtained using the self.get_all_label_spaces():
            {'mat': 1, 'time': 1}
            {'mat': 0, 'time': 1}
            {'mat': 1, 'time': 2}
            {'mat': 0, 'time': 2}
        To get the plotted result at the time_step number 2, use: self.plot_contour("time", [1])
        """
        self._evaluate_result()
        if not sys.warnoptions:
            import warnings
            warnings.simplefilter("ignore")
        plotter = pv.Plotter()
        mesh = self._evaluator._model.metadata.meshed_region
        grid = mesh.grid
        nan_color = "grey"
        rescoper = _Rescoper(mesh, self.result_fields_container[0].location, 
                             self.result_fields_container[0].component_count) #location will be the same on all fields
        if (self.result_fields_container.__len__() == 1):
            field = rescoper.rescope(self.result_fields_container[0])
            plotter.add_mesh(grid, scalars = field, opacity=1.0, nan_color=nan_color, 
                              stitle = self.result_fields_container[0].name, show_edges=True)
            
        else:
            label_spaces = []
            for opt_id in option_id:
                i = 0
                while i < self.result_fields_container.__len__():
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
            for label in label_spaces:
                field_to_rescope = self.result_fields_container._get_entries(label)
                if field_to_rescope is None:
                    raise Exception("The label " + label.__str__() + " does not exist in the fields container.")
                name = self.result_fields_container[0].name.split("_")[0]
                field = rescoper.rescope(field_to_rescope)
                plotter.add_mesh(grid, scalars = field, nan_color=nan_color, stitle = name, show_edges=True)
        plotter.add_axes()
        plotter.show()
        
    
        
    def plot_chart(self):
        """Plot the minimum/maximum result values over time 
        if the time_freq_support contains several time_steps 
        (for example: transient analysis)"""
        tfq = self._model.metadata.time_freq_support
        timeids = list(range(1,tfq.n_sets+1))
        self._result_operator.inputs.time_scoping.connect(timeids)
        self._evaluate_result_forced()
        time_field = tfq.frequencies
        normOp = Operator("norm_fc")
        minmaxOp = Operator("min_max_fc")
        normOp.inputs.fields_container.connect(self.result_fields_container)
        minmaxOp.inputs.connect(normOp.outputs)
        fieldMin = minmaxOp.outputs.field_min()
        fieldMax = minmaxOp.outputs.field_max()
        pyplot.plot(time_field.data,fieldMax.data,'r',label='Maximum')
        pyplot.plot(time_field.data,fieldMin.data,'b',label='Minimum')
        pyplot.xlabel("time (s)")
        pyplot.ylabel(self._result_operator.name + fieldMin.unit)
        substr = self.result_fields_container[0].name.split("_")
        pyplot.title( substr[0] + ": min/max values over time")
        pyplot.legend()
        
        
    def is_complex_result(self):
        """Returns True if the result contains complex frequencies."""
        tfq = self.result_fields_container.time_freq_support
        if (tfq.complex_frequencies != None):
            return True
        else:
            return False
        
        
        
    