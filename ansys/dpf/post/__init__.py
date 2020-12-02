from ansys.dpf.post.common import Grouping as grouping
from ansys.dpf.core.common import locations
from ansys.dpf.post.common import _AvailableKeywords as available_keywords
from ansys.dpf.post.post_utility import solution, build_docs
from ansys import dpf

"""Post-processing module. Using Data Processing Framework.
Allow to create a result object, then use it to get wanted results.

Example
-----
from ansys.dpf import post
result = post.result("file.rst")
disp = result.nodal_displacement() 

"""

dpf.core.start_local_server()

    



