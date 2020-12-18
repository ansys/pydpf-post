"""Module containing the common tools for a better usage of the DPF-Post module."""

from enum import Enum


# class ElShapes(Enum):
#     """Class with Enum inheritance. Must be used to 
#     describe the element shape when API allows it.
    
#     Example
#     -----
#     from ansys.dpf import post
#     result = post.result("file.rst")
#     disp = result.elemental_stress(element_shape = post.el_shape.shell)
#     """
#     solid = 1
#     beam = 2
#     shell = 3
#     shell_top = 4
#     shellmid = 5
#     shell_bottom = 6
    
    
class Grouping():
    """Class with Enum inheritance. Must be used to 
    describe a grouping command when the API allows it.
    
    Examples
    --------
    >>> from ansys.dpf import post
    >>> solution = post.solution("file.rst")
    >>> disp = solution.elemental_stress(element_shape = post.grouping.by_el_shape)
    """
    by_el_shape = "elshape"
    by_material = "mat"
    by_body = "body"
    
    
class _AvailableKeywords():
    """Contains all the keywords that can be used inside of 
    a method from a post.solution(file_path) object. 
    
    In order to view the complete list of available keywords, use:
        post.print_available_keywords()
    """
    location = "location"
    node_scoping = "node_scoping"
    element_scoping = "element_scoping"
    named_selection = "named_selection"
    time = "time"
    set = "set"
    grouping = "grouping"
    _phase = "phase"
    _subresult = "subresult"
    mapdl_grouping = "mapdl_grouping"
    time_scoping = "time_scoping"
    
    def __str__(self):
        txt = ""
        for attr in dir(_AvailableKeywords):
            if not attr.startswith("__") and not attr.startswith("_") and not callable(getattr(_AvailableKeywords, attr)):
                txt += attr
                txt += ": "
                txt += self._description_mapping(attr)
                txt += "\n"
        return txt
    
    def _description_mapping(self, attr_name):
        if (attr_name == self.location):
            return "str. Use post.locations.(...) as helper."
        if (attr_name == self.node_scoping):
            return "list, int or dpf.core.Scoping"
        if (attr_name == self.element_scoping):
            return "list, int or dpf.core.Scoping"
        if (attr_name == self.time_scoping):
            return "list, int or dpf.core.Scoping"
        if (attr_name == self.named_selection):
            return "str. Name of named_selection."
        if (attr_name == self.time):
            return "float"
        if (attr_name == self.set):
            return "int"
        if (attr_name == self.mapdl_grouping):
            return "int. Write 186 to get mapdl_elements solid_186."
        if (attr_name == self.grouping):
            return "str. Use post.grouping.(...) as helper."
    
    
class _AnalysisType():
    """Contains Python analysis type names. For developers usage."""
    static = "static"
    modal = "modal"
    harmonic = "harmonic"
    transient = "transient"
    
    
class _PhysicsType():
    """Contains Python physics type names. For developers usage."""
    mecanic = "mecanic"
    thermal = "thermal"
    

        
