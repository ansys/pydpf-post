"""Module containing the ResultEvaluation class that
will compute the information collected in the ResultData 
class (which is a fields container wrapper).

This class will only be dedicated to the fields container 
computation."""

import numpy as np
from collections import OrderedDict
from ansys import dpf
from ansys.dpf.core import Operator
from ansys.dpf.core.common import types, locations
from ansys.dpf.post.common import Grouping, _AvailableKeywords


class ResultEvaluator:
    """This object will make the evaluation of the fields container 
    wrapped in the ResultData object."""
    def __init__(self, operator_name: str, data_sources, model,
                 elem_average, 
                 location: str = None, element_scoping = None, 
                 node_scoping = None, named_selection = None, 
                 time = None, grouping = None, phase = None, subresult = None, 
                 mapdl_grouping = None, set = None, time_scoping = None):
        
        self._model = model
        self._chained_operators = OrderedDict() #dictionary containing (key = Operator.name, [Operator, description])
        self.subresult = subresult
        
        if (subresult != None):
            operator_name += subresult
        self._result_operator = dpf.core.Operator(operator_name)
        self._result_operator.inputs.connect(data_sources)
        
        if (location != None):
            self._result_operator.inputs.requested_location.connect(location)
        if (element_scoping != None and node_scoping != None):
            raise Exception("Impossible to use both element_scoping and node_scoping.")
        if (set != None and time != None) or (set != None and time_scoping != None) or (time != None and time_scoping != None):
            raise Exception("Set, time and time_scoping keyword arguments must be used independently.")
        if (element_scoping != None):
            scoping = self._compute_scoping(element_scoping, locations.elemental)
            self._result_operator.inputs.mesh_scoping.connect(scoping)
        if (node_scoping != None):
            scoping = self._compute_scoping(node_scoping, locations.nodal)
            self._result_operator.inputs.mesh_scoping.connect(scoping)
        if (time_scoping != None):
            t_scoping = self._compute_scoping(time_scoping)
            self._result_operator.inputs.time_scoping.connect(t_scoping)  
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
            mesh_provider.inputs.data_sources.connect(data_sources)
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
            scop_by_prop_op.inputs.data_sources.connect(data_sources)
            scop_by_prop_op.inputs.property_id.connect(mapdl_grouping)
            scop_by_prop_op.inputs.requested_location.connect(locations.elemental)
            self._result_operator.inputs.mesh_scoping.connect(scop_by_prop_op.outputs.mesh_scoping)
            self._chained_operators[scop_by_prop_op.name] = """This operator will compute a scoping from a mapdl elemen type id. Its output (mesh_scoping) will be connected with the mesh_scoping input of the result operator."""
        if (set != None):
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
        # outside post-processing instruction         
        if elem_average:
            self._elemental_nodal_to_elemental_result()
                
    def _compute_scoping(self, in_scoping, in_location = None):
        out_scoping = in_scoping
        if not isinstance(in_scoping, dpf.core.scoping.Scoping):
            out_scoping = dpf.core.Scoping()
            if in_location != None:
                out_scoping.location = in_location
            if isinstance(in_scoping, list):
                out_scoping.ids = in_scoping
            elif isinstance(in_scoping, int):
                l = [in_scoping]
                out_scoping.ids = l
            else:
                error_scoping = "Only dpf.core.Scoping list or int are supported as scoping."
                raise Exception(error_scoping)
        return out_scoping
            
    def _check_if_scoping(self, node_scoping, element_scoping):
        if (node_scoping != None) or (element_scoping != None):
            txt = "Keywords " + _AvailableKeywords.element_scoping + "/" + _AvailableKeywords.node_scoping + " can not be used with " + _AvailableKeywords.grouping + "/" + _AvailableKeywords.named_selection + " ones."
            raise Exception(txt)
      
    def _check_if_several_grouping(self, grouping, mapdl_grouping):
        if (grouping != None) and (mapdl_grouping != None):
            raise Exception("Both keywords grouping and mapdl_grouping can not be used simultaneously.")
           
    def _get_evaluation_with_sweeping_phase(self, phase):
        """Connects needed operator to compute the result regarding the specified phase."""
        sweeping_phase_op = dpf.core.Operator("sweeping_phase_fc")
        sweeping_phase_op.inputs.fields_container.connect(self._result_operator.outputs.fields_container)
        sweeping_phase_op.inputs.angle.connect(phase)
        sweeping_phase_op.inputs.unit_name.connect("deg")
        self._chained_operators[sweeping_phase_op.name] = """This operator will compute the result at a given phase (when result has complex values)."""
        self._result_operator = sweeping_phase_op      
        
    def _elemental_nodal_to_elemental_result(self):
        avg = Operator("to_elemental_fc")
        fc_to_connect = self._result_operator.outputs.fields_container
        avg.inputs.fields_container.connect(fc_to_connect)
        self._result_operator = avg
            
    def evaluate_result(self):
        """Re-evaluation of the result."""
        result_fc = self._result_operator.get_output(0, types.fields_container)
        return result_fc