"""This module contains the super class of the 
stress/strain/temperature/displacement objects."""

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.core.common import locations
from ansys.dpf.post.result_data import ResultData

class Result:
    def __init__(self, data_sources, model, **kwargs):
        self._data_sources = data_sources
        self._model = model
        
        self._location = locations.nodal
        self._element_scoping = None
        self._node_scoping = None
        self._named_selection = None
        self._grouping = None
        self._mapdl_grouping = None
        self._time_scoping = None
        if _AvailableKeywords.location in kwargs:
            self._location = kwargs.pop(_AvailableKeywords.location)
        if _AvailableKeywords.element_scoping in kwargs:
            self._element_scoping = kwargs.pop(_AvailableKeywords.element_scoping)
        if _AvailableKeywords.node_scoping in kwargs:
            self._node_scoping = kwargs.pop(_AvailableKeywords.node_scoping)
        if _AvailableKeywords.named_selection in kwargs:
            self._named_selection = kwargs.pop(_AvailableKeywords.named_selection)
        if _AvailableKeywords.grouping in kwargs:
            self._grouping = kwargs.pop(_AvailableKeywords.grouping)
        if _AvailableKeywords.mapdl_grouping in kwargs:
            self._mapdl_grouping = kwargs.pop(_AvailableKeywords.mapdl_grouping)
        if _AvailableKeywords.time_scoping in kwargs:
            self._time_scoping = kwargs.pop(_AvailableKeywords.time_scoping)
            
        if (kwargs.__len__() > 0):
            txt = "Following keyword arguments could not be computed: \n"
            for key in kwargs:
                txt += "- " + key + "\n"
            txt += "\n"
            txt += "Use 'print(post.available_keywords())' to see which keywords can be used, \n"
            txt += "or try to use the keyword while calling the subresult (if stress = solution.stress(), \n"
            txt += "try stress_x = solution.stress(MY_KEYWORDS=KEYWORDS_VALUE))."
            raise Exception(txt)
            
    def __str__(self):
        txt = "Object properties are: \n"
        for attr, value in vars(self).items():
            if not (attr.startswith('__') or attr.startswith('_data_sources') or attr.startswith('_model')):
                if value is not None:
                    if attr.startswith('_'):
                        attr = attr[1:]
                    if isinstance(value, list):
                        value_temp = "["
                        for val in value:
                             value_temp += str(val) + ","
                        value_temp = value_temp[:value_temp.__len__()-1]
                        value_temp += "]"
                        value = value_temp
                    txt += " - "+ attr
                    txt += " = " + value + "\n"
        return txt
    
    @property
    def location(self):
        """Location property. Can be set. 
        Type: str. 
        post.locations enum can be used.
        
        Example
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.location = post.locations.elemental
        """
        return self._location

    @location.setter
    def location(self, value):
        self._location = value
        
    @property
    def element_scoping(self):
        """Elemental scoping property. Can be set. 
        Type: list, or dpf.core.Scoping. 
        
        Example
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.element_scoping = [1, 4]
        """
        return self._element_scoping

    @element_scoping.setter
    def element_scoping(self, value):
        self._element_scoping = value
        
    @property
    def node_scoping(self):
        """Nodal scoping property. Can be set. 
        Type: list, or dpf.core.Scoping. 
        
        Example
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.node_scoping = [1, 4]
        """
        return self._node_scoping

    @node_scoping.setter
    def node_scoping(self, value):
        self._node_scoping = value
        
    @property
    def named_selection(self):
        """Named selection property. Can be set. 
        Type: str.
        The named selection name must be written following upper case.
        
        Example
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.node_scoping = 'SELECTION'
        """
        return self._named_selection

    @named_selection.setter
    def named_selection(self, value):
        self._named_selection = value
        
    @property
    def grouping(self):
        """Grouping property. Can be set. 
        Type: str.
        post.grouping enum can be used.
        
        Example
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.node_scoping = post.grouping.by_el_shape
        """
        return self._grouping

    @grouping.setter
    def grouping(self, value):
        self._grouping = value
        
    @property
    def time_scoping(self):
        """Time scoping property. Can be set. 
        Type: list, or dpf.core.Scoping. 
        
        Example
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.time_scoping = [1, 4]
        """
        return self._time_scoping

    @time_scoping.setter
    def time_scoping(self, value):
        self._time_scoping = value
        
    @property
    def mapdl_grouping(self):
        """Mapdl grouping property. Can be set. 
        Type: int.
        
        Example to get only solid 186 mapdl element type:
        -----
        from ansys.dpf import post
        solution = post.solution("file.rst")
        stress = solution.stress()
        stress.mapdl_grouping = 186
        """
        return self._mapdl_grouping

    @mapdl_grouping.setter
    def mapdl_grouping(self, value):
        self._mapdl_grouping = value
        
        
            
    def _get_result_data(self, operator_name, data_sources, model, **kwargs):
        """This method checks the keyword arguments that are 
        specified while calling a subresult method.
        
        The arguments that can be set at this point are:
            - time
            - set
            - phase (if complex result)
        """
        #write correct arguments regarding location
        b_elem_average = False
        location_to_compute = self._location
        if self._location == locations.elemental_nodal:
            location_to_compute = locations.elemental
        if self._location == locations.elemental:
            b_elem_average = True
            
        #user keywords
        time = None
        phase = None
        set = None
        if _AvailableKeywords.phase in kwargs:
            phase = kwargs[_AvailableKeywords.phase]    
        if _AvailableKeywords.time in kwargs:
            time = kwargs[_AvailableKeywords.time]
        if _AvailableKeywords.set in kwargs:
            set = kwargs[_AvailableKeywords.set]
            
        #!TODO check the previous set keywords ?
            
        #program keywords
        subresult=None
        if _AvailableKeywords.subresult in kwargs:
            subresult = kwargs[_AvailableKeywords.subresult]
        return ResultData(operator_name=operator_name, data_sources=data_sources,
                          model=model, elem_average=b_elem_average, location=location_to_compute, 
                          element_scoping=self._element_scoping, node_scoping=self._node_scoping, 
                          named_selection=self._named_selection,
                          time=time, grouping=self._grouping, phase=phase, 
                          subresult=subresult, mapdl_grouping=self._mapdl_grouping, set=set, 
                          time_scoping=self._time_scoping)
