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

from ansys.dpf.post.result import result
from ansys.dpf.post.common import ElShapes as el_shape
from ansys.dpf.post.common import Grouping as grouping
from ansys import dpf

"""Post-processing module. Using Data Processing Framework.
Allow to create a result object, then use it to get wanted results.

Example
-----
from ansys import post
result = post.result("file.rst")
disp = result.nodal_displacement() 

"""

dpf.core.start_local_server()