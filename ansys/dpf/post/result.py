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

from ansys.dpf.core.model import Model
from ansys.dpf.post.flat_result import StaticAnalysisResult, ModalAnalysisResult, HarmonicAnalysisResult, TransientAnalysisResult

def result(data_sources=None, channel=None):
    """Return a Result object which can provide information on a given set, on a given scoping...
    
    Parameters
    ----------
    Can be a filepath to the file you want to open, or a dpf.DataSources().
        
    Examples
    --------
    result = post.Result("file.rst")
    """
    
    _model = Model(data_sources, channel)
    # self._base = BaseService(self._channel)
    # self.metadata = Metadata(data_sources, channel)
    # self.results = Results(self)
    
    
    analysis_type = _model.metadata.result_info.analysis_type
    data_sources = _model.metadata.data_sources
    if(analysis_type == 'static'):
        return StaticAnalysisResult(data_sources, _model.metadata)
    elif (analysis_type == 'modal'):
        return ModalAnalysisResult(data_sources, _model.metadata)
    elif (analysis_type == 'harmonic'):
        return HarmonicAnalysisResult(data_sources, _model.metadata)
    elif (analysis_type == 'transient'):
        return TransientAnalysisResult(data_sources, _model.metadata)
    else:
        raise Exception("Unknown model.metadata.result_info.")
        
