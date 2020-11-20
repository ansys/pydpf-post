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

from enum import Enum

class ElShapes(Enum):
    solid = 1
    beam = 2
    shell = 3
    shell_top = 4
    shellmid = 5
    shell_bottom = 6
    
class Grouping(Enum):
    by_el_shape = 1
    by_el_type = 2
    by_material = 3
    by_body = 4


def _map_property_name(property_enum_value):
    if (property_enum_value == 1):
        return "elshape"
    elif (property_enum_value == 2):
        return "eltype"
    elif (property_enum_value == 3):
        return "mat"
    elif (property_enum_value == 4):
        return "body"
        
