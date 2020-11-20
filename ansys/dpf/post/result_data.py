##########################################################################
#                                                                        #
#          Copyright (C) 2020 ANSYS Inc.  All Rights Reserved            #
#                                                                        #
# This file contains proprietary software licensed from ANSYS Inc.       #
# This header must remain in any source code despite modifications or    #
# enhancements by any party.                                             #
#                                                                        #
##########################################################################
# Version: 1.0                                                           #
# Author(s): C.Bellot/R.Lagha/L.Paradis                                  #
# contact(s): ramdane.lagha@ansys.com                                    #
##########################################################################

"""Module containing the _ResultData class ("result object") that
user will be able to compute through the DPF Flat API."""

import numpy as np
from ansys import dpf
from ansys.dpf.core import Operator
from ansys.dpf.core.common import types
from ansys.dpf.post.common import ElShapes, Grouping

class ResultData:
    def __init__(self, operator_name: str, data_sources, 
                 location: str = None, element_scoping = None, 
                 node_scoping = None, named_selection = None, 
                 el_shape = None, time_step = None, 
                 grouping = None, phase = None):
        self.result_fields_container = None
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
                scoping.location = 'Elemental'
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
                scoping.location = 'Nodal'
                if isinstance(node_scoping, list):
                    scoping.ids = node_scoping
                elif isinstance(node_scoping, np.array):
                    scoping.ids = node_scoping.tolist()
                else:
                    raise Exception("Only dpf.core.Scoping, numpy.array or list are supported as scoping.")
            self._result_operator.inputs.mesh_scoping.connect(scoping)
        if (named_selection != None):
            ns_op = dpf.core.Operator("scoping_provider_by_ns")
            ns_op.inputs.data_sources.connect(data_sources)
            ns_op.inputs.named_selection_name.connect(named_selection)
            self._result_operator.inputs.mesh_scoping.connect(ns_op.outputs.mesh_scoping)
        # #!TODO time_step and substep, el_shape, grouping, set 
        # if (grouping != None):
        #     scop_op = Operator("scoping::by_property")
        # get a mapping between enums from dataProcessingCore.dll and python. Grpc?
        #     scop_op.inputs.property_name.connect()
        #     self._result_operator.inputs.mesh_scoping.connect(scop_op.outputs)
        if (phase != None):
            if not (isinstance(phase, float) or isinstance(phase, int)):
                raise Exception("Phase angle must be an int or a float number. The angle unit is degree.")
            self._get_evaluation_with_sweeping_phase(phase)
        
    def _evaluate_result(self):
        if (self.result_fields_container == None):
            result_fc = self._result_operator.get_output(0, types.fields_container)
            self.result_fields_container = result_fc
        
    def _get_evaluation_with_sweeping_phase(self, phase):
        sweeping_phase_op = dpf.core.Operator("sweeping_phase_fc")
        sweeping_phase_op.inputs.fields_container.connect(self._result_operator.outputs.fields_container)
        sweeping_phase_op.inputs.angle.connect(phase)
        sweeping_phase_op.inputs.unit_name.connect("deg")
        self._result_operator = sweeping_phase_op
        
        
    def num_fields(self):
        self._evaluate_result()
        return self.result_fields_container.__len__()
    
    def data_at_field(self, field_index: int = 0):
        self._evaluate_result()
        return self.result_fields_container[field_index].data
    
    def __getitem__(self, field_index: int = 0):
        self._evaluate_result()
        return self.result_fields_container[field_index]
    
    #!TODO
    def scoping_at_field(self, field_index: int = 0):
        self._evaluate_result()
        raise Exception("Not implemented yet.")
        
    def _min_max(self, pin):
        self._evaluate_result()
        max_operator = Operator("min_max_fc")
        max_operator.inputs.connect(self.result_fields_container)
        result = max_operator.get_output(pin, types.field)
        return result
        
    def max(self):
        return self._min_max(1)
    
    def max_data(self):
        return self._min_max(1).data
    
    def max_data_at_field(self, field_index: int = 0):
        return self._min_max(1).data[field_index]
    
    def min(self):
        return self._min_max(0)
    
    def min_data(self):
        return self._min_max(0).data
    
    def min_data_at_field(self, field_index: int = 0):
        return self._min_max(0).data[field_index]
    
    #!TODO
    def eplot(self):
        self._evaluate_result()
        raise Exception("Not implemented yet.")
        
        
    