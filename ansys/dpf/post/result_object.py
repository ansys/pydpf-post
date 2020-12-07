"""This module contains the super class of the 
stress/strain/temperature/displacement objects."""

from ansys.dpf.core.common import locations
from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.result_definition import Definition

class Result:
    def __init__(self, data_sources, model, **kwargs):
        self._data_sources = data_sources
        self._model = model
        
        self.definition = Definition(**kwargs)
        
        
    def __str__(self):
        txt = self.definition.__str__()
        return txt
    
            
    def _get_result_data(self, operator_name, data_sources, model, subresult=None):
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
                          model=model, elem_average=b_elem_average, location=location_to_compute, 
                          element_scoping=self.definition.element_scoping, node_scoping=self.definition.node_scoping, 
                          named_selection=self.definition.named_selection,
                          time=self.definition.time, grouping=self.definition.grouping, 
                          subresult=subresult, mapdl_grouping=self.definition.mapdl_grouping, set=self.definition.set, 
                          time_scoping=self.definition.time_scoping)
