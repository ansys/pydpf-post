"""Module containing the _ResultData class ("result object") that
user will be able to compute through the DPF Post API."""

import numpy as np
import pyvista as pv
import matplotlib.pyplot as pyplot
import os
from collections import OrderedDict
from ansys import dpf
from ansys.dpf.core import Operator
from ansys.dpf.core.common import types, locations
from ansys.dpf.post.common import ElShapes, Grouping, _AvailableKeywords


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
        - el_shape
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
    
    def __init__(self, operator_name: str, data_sources, model, instance, 
                 location: str = None, element_scoping = None, 
                 node_scoping = None, named_selection = None, 
                 el_shape = None, time = None, 
                 grouping = None, phase = None, subresult = None, 
                 mapdl_grouping = None, set = None):
        self._model = model
        self.result_fields_container = None
        self._chained_operators = OrderedDict() #dictionary containing (key = Operator.name, [Operator, description])
        
        self._data_sources = data_sources
        if (subresult != None):
            operator_name += subresult
        self._result_operator = dpf.core.Operator(operator_name)
        self._result_operator.inputs.connect(data_sources)

        if (location != None):
            self._result_operator.inputs.requested_location.connect(location)
        if (element_scoping != None and node_scoping != None):
            raise Exception("Impossible to use both element_scoping and node_scoping.")
        if (element_scoping != None):
            scoping = element_scoping
            if not isinstance(element_scoping, dpf.core.scoping.Scoping):
                scoping = dpf.core.Scoping()
                scoping.location = locations.elemental
                if isinstance(element_scoping, list):
                    scoping.ids = element_scoping
                elif isinstance(element_scoping, np.array):
                    scoping.ids = element_scoping.tolist()
                else:
                    raise Exception("Only dpf.core.Scoping, numpy.array or list are supported as scoping.")
            self._result_operator.inputs.mesh_scoping.connect(scoping)
        if (node_scoping != None):
            scoping = node_scoping
            if not isinstance(node_scoping, dpf.core.scoping.Scoping):
                scoping = dpf.core.Scoping()
                scoping.location = locations.nodal
                if isinstance(node_scoping, list):
                    scoping.ids = node_scoping
                elif isinstance(node_scoping, np.array):
                    scoping.ids = node_scoping.tolist()
                else:
                    raise Exception("Only dpf.core.Scoping, numpy.array or list are supported as scoping.")
            self._result_operator.inputs.mesh_scoping.connect(scoping)
        #chained before the result_operator
        if (named_selection != None):
            self._check_if_scoping(node_scoping, element_scoping)
            ns_op = Operator("scoping_provider_by_ns")
            ns_op.inputs.data_sources.connect(data_sources)
            ns_op.inputs.named_selection_name.connect(named_selection)
            self._result_operator.inputs.mesh_scoping.connect(ns_op.outputs.mesh_scoping)
            self._chained_operators[ns_op.name] = """This operator will compute a scoping from a named selection name. Its output (mesh_scoping) will be connected with the mesh_scoping input of the result operator."""
        if (grouping != None):
            self._check_if_scoping(node_scoping, element_scoping)
            # part for grouping by_material/by_el_shape
            mesh_provider = Operator("MeshProvider")
            mesh_provider.inputs.data_sources.connect(self._data_sources)
            scop_op = Operator("scoping::by_property")
            scop_op.inputs.mesh.connect(mesh_provider.outputs.mesh) #default output location is elemental
            if (location != None):
                if(location == locations.nodal):
                    scop_op.inputs.requested_location.connect(locations.nodal)
                else:
                    scop_op.inputs.requested_location.connect(locations.elemental)
            else:
                scop_op.inputs.requested_location.connect(locations.nodal)
            if (grouping == Grouping.by_material) or (grouping == Grouping.by_body):
                scop_op.inputs.label1.connect("mat")
            elif(grouping == Grouping.by_el_shape):
                 scop_op.inputs.label1.connect("elshape")            
            else:
                raise Exception("Grouping impossible. Use the keyword argument as: grouping = grouping.by_el_shape, grouping = grouping.by_material...")
            self._result_operator.inputs.mesh_scoping.connect(scop_op.outputs.mesh_scoping)
            self._chained_operators[scop_op.name] = """This operator will compute a scoping from a grouping option. Its output (mesh_scoping) will be connected with the mesh_scoping input of the result operator."""
        if (mapdl_grouping != None):
            self._check_if_scoping(node_scoping, element_scoping)
            # part for grouping by_el_type
            #!TODO, replace grouping.by_el_type_mapdl... by grouping.by_el_type (general case)
            scop_by_prop_op = Operator("scoping_provider_by_prop")
            scop_by_prop_op.inputs.property_name.connect("mapdl_element_type")
            scop_by_prop_op.inputs.data_sources.connect(self._data_sources)
            scop_by_prop_op.inputs.property_id.connect(mapdl_grouping)
            scop_by_prop_op.inputs.requested_location.connect(locations.elemental)
            self._result_operator.inputs.mesh_scoping.connect(scop_by_prop_op.outputs.mesh_scoping)
            self._chained_operators[scop_by_prop_op.name] = """This operator will compute a scoping from a mapdl elemen type id. Its output (mesh_scoping) will be connected with the mesh_scoping input of the result operator."""
        if (set != None):
            if (time != None):
                raise Exception("'time' and 'set' keyword can not be used simultaneously.")
            if not isinstance(set, int):
                raise Exception("Set argument must be an int value.")
            time_scoping = dpf.core.Scoping()
            time_scoping.ids = [set]
            self._result_operator.inputs.time_scoping.connect(time_scoping)
        #add the result operator
        self._chained_operators[self._result_operator.name] = "Result operator. Compute the wanted result"
        #chained after the result_operator
        if (phase != None):
            self._get_evaluation_with_sweeping_phase(phase)
        if (time != None):
            if (set != None):
                raise Exception("'time' and 'set' keyword can not be used simultaneously.")
            if not isinstance(time, float):
                raise Exception("Time argument must be a float value.")
            time_scoping = dpf.core.Scoping()
            tfq = self._model.metadata.time_freq_support
            data = tfq.frequencies.data
            temp_array = np.array([])
            for d in data:
                temp_array = np.append(temp_array, round(d, 5))
            if time in temp_array:
                index = tfq.get_cumulative_index(freq = time)
                time_scoping.ids = [index + 1]
                self._result_operator.inputs.time_scoping.connect(time_scoping)
            else:
                # centroid when time value is between to time steps
                lower_index = tfq.get_cumulative_index(freq = time)
                time_scoping.ids = [lower_index+1, lower_index+2]
                self._result_operator.inputs.time_scoping.connect(time_scoping)
                centroid_op = Operator("centroid")
                time1 = tfq.get_frequency(cumulative_index=lower_index)
                time2 = tfq.get_frequency(cumulative_index=(lower_index+1))
                factor = (time - time1) / (time2 - time1)
                centroid_op.inputs.factor.connect(factor)
                outp = self._result_operator.outputs.fields_container()
                fieldA = outp[0]
                fieldB = outp[1]
                centroid_op.inputs.fieldA.connect(fieldA)
                centroid_op.inputs.fieldB.connect(fieldB)
                self._chained_operators[centroid_op.name] = "This operator will compute the centroid of two fields obtained with a time scoping containing two times."
                forward_op = Operator("forward_fc")
                forward_op.inputs.fields.connect(centroid_op.outputs.field)
                self._result_operator = forward_op
                
            
    def __str__(self):
        self._evaluate_result()
        name = self.result_fields_container[0].name.split("_")
        txt = "%s result.\n\n" % name[0].capitalize()
        txt += "The result is computed thanks to dpf.core.Operator objects, "
        txt += "that are chained together regarding the following list: \n"
        for key, val in self._chained_operators.items():
            txt += "- %s: " % key
            txt += val
            txt += "\n"
        return txt
    
    def _check_if_scoping(self, node_scoping, element_scoping):
        if (node_scoping != None) or (element_scoping != None):
            txt = "Keywords " + _AvailableKeywords.element_scoping + "/" + _AvailableKeywords.node_scoping + " can not be used with " + _AvailableKeywords.grouping + "/" + _AvailableKeywords.named_selection + " ones."
            raise Exception(txt)
      
    def _check_if_several_grouping(self, grouping, mapdl_grouping):
        if (grouping != None) and (mapdl_grouping != None):
            raise Exception("Both keywords grouping and mapdl_grouping can not be used simultaneously.")
        
    def _evaluate_result(self):
        """First evaluation of the result."""
        if (self.result_fields_container == None):
            result_fc = self._result_operator.get_output(0, types.fields_container)
            self.result_fields_container = result_fc
            
    def _evaluate_result_forced(self):
        """Re-evaluation of the result."""
        result_fc = self._result_operator.get_output(0, types.fields_container)
        self.result_fields_container = result_fc
        
    def _get_evaluation_with_sweeping_phase(self, phase):
        """Connects needed operator to compute the result regarding the specified phase."""
        sweeping_phase_op = dpf.core.Operator("sweeping_phase_fc")
        sweeping_phase_op.inputs.fields_container.connect(self._result_operator.outputs.fields_container)
        sweeping_phase_op.inputs.angle.connect(phase)
        sweeping_phase_op.inputs.unit_name.connect("deg")
        self._chained_operators[sweeping_phase_op.name] = """This operator will compute the result at a given phase (when result has complex values)."""
        self._result_operator = sweeping_phase_op
        
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
    
    #!TODO
    def scoping_at_field(self, field_index: int = 0):
        """Returns the scoping of the result."""
        self._evaluate_result()
        raise Exception("Not implemented yet.")
        
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
    
    
    def plot_contour(self):
        """Plot the contour result on its mesh support. The obtained figure depends on the 
        support (can be a meshed_region or a time_freq_support).
        If transient analysis, plot the last result."""
        self._evaluate_result()
        plotter = pv.Plotter()
        mesh_provider = Operator("MeshProvider")
        mesh_provider.inputs.data_sources.connect(self._data_sources)
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
        
    
    def _plot_contour_try(self):
        """Plot contour, does not work yet (plot_contour() must be used instead)."""
        self._evaluate_result()
        plotter = pv.Plotter()
        #comment savoir si je dois plotter en fct du maillage ou du tfq ? 
        #quid des fc qui contiennent des fields séparés pour un même modèle ? (par ex: fct de el_shape)
        #-> Operator("GetSupportFromField") will return the mesh of a given field
        #-> plotter.add_mesh(scalars...)
        #quid si scoping ? Pour l'instant pas ok
        if (self.result_fields_container.__len__() == 1):
            grid = self._model.metadata.meshed_region.grid
            field = self.result_fields_container[0]
            plotter.add_mesh(grid, scalars = field.data, stitle = field.name, show_edges=True)
        else:
            name = self.result_fields_container[0].name
            for field in self.result_fields_container:
                if (field.name == name):
                    mesh_op = Operator("GetSupportFromField")
                    mesh_op.inputs.field.connect(field)
                    mesh = mesh_op.outputs.mesh()
                    grid = mesh.grid
                    plotter.add_mesh(grid, scalars = field.data, stitle = field.name, show_edges=True)
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
        
        
        
    