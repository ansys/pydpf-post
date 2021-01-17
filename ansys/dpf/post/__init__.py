from ansys.dpf.post.common import Grouping as grouping
from ansys.dpf.core.common import locations
from ansys.dpf.post.post_utility import load_solution, print_available_keywords
from ansys import dpf

from ansys.dpf.post.misc import Report

"""Post-processing module. Using Data Processing Framework.
Allow to create a result object, then use it to get wanted results.

Examples
--------
>>> from ansys.dpf import post
>>> solution = post.solution("file.rst")
>>> disp = solution.nodal_displacement() 

"""

#dpf.core.start_local_server()

    
