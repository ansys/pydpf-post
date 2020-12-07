"""Module containing the methode to instantiate the result object. Initialization of post objects."""

from ansys.dpf.core.model import Model 
 
from ansys.dpf.post.common import _AnalysisType, _AvailableKeywords
from ansys.dpf.post.static_analysis import StaticAnalysisSolution
from ansys.dpf.post.modal_analysis import ModalAnalysisSolution
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.transient_analysis import TransientAnalysisSolution
        

def load_solution(data_sources):
    """Return a Result object which can provide information on a given set, on a given scoping...
    
    Parameters
    ----------
    str
        Can be a filepath to the file you want to open, or a dpf.core.DataSources().
        
    Examples
    --------
    solution = post.solution("file.rst")
    """
    _model = Model(data_sources)
    analysis_type = _model.metadata.result_info.analysis_type
    data_sources = _model.metadata.data_sources
    if(analysis_type == _AnalysisType.static):
        return StaticAnalysisSolution(data_sources, _model)
    elif (analysis_type == _AnalysisType.modal):
        return ModalAnalysisSolution(data_sources, _model)
    elif (analysis_type == _AnalysisType.harmonic):
        return HarmonicAnalysisSolution(data_sources, _model)
    elif (analysis_type == _AnalysisType.transient):
        return TransientAnalysisSolution(data_sources, _model)
    else:
        raise Exception("Unknown model.metadata.result_info.")
        
        
def print_available_keywords():
    txt = _AvailableKeywords().__str__()
    print(txt)
     