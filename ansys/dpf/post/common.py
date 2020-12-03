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
    
    
class Grouping(Enum):
    """Class with Enum inheritance. Must be used to 
    describe a grouping command when the API allows it.
    
    Examples
    --------
    >>> from ansys.dpf import post
    >>> solution = post.solution("file.rst")
    >>> disp = solution.elemental_stress(element_shape = post.grouping.by_el_shape)
    """
    by_el_shape = 1
    by_material = 2
    by_body = 3
    
    
class _AvailableKeywords():
    """Contains all the keywords that can be used inside of 
    a method from a post.solution(file_path) object. 
    
    In order to view the complete list of available keywords, use:
        print(post.available_kewords())
    """
    location = "location"
    node_scoping = "node_scoping"
    element_scoping = "element_scoping"
    named_selection = "named_selection"
    time = "time"
    # substep = "substep"
    set = "set"
    grouping = "grouping"
    phase = "phase"
    subresult = "subresult"
    mapdl_grouping = "mapdl_grouping"
    time_scoping = "time_scoping"
    
    def __str__(self):
        txt = ""
        for attr in dir(_AvailableKeywords):
            if not attr.startswith("__") and not callable(getattr(_AvailableKeywords, attr)):
                txt += attr
                txt += "\n"
        return txt
    
    
class _AnalysisType():
    """Contains Python analysis type names. For developers usage.
    """
    static = "static"
    modal = "modal"
    harmonic = "harmonic"
    transient = "transient"
    

def _map_property_name(property_enum_value):
    """Map to get property name from an enum value. For developers usage.
    """
    if (property_enum_value == 1):
        return "elshape"
    elif (property_enum_value == 2):
        return "eltype"
    elif (property_enum_value == 3):
        return "mat"
    elif (property_enum_value == 4):
        return "body"
        
