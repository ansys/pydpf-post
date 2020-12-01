"""Module containing the methode to instantiate the result object. Initialization of post objects."""

from ansys.dpf.core.model import Model 
import os
from ansys import dpf
 
from ansys.dpf.post.common import _AnalysisType
from ansys.dpf.post.static_analysis import StaticAnalysisResult
from ansys.dpf.post.modal_analysis import ModalAnalysisResult
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisResult
from ansys.dpf.post.transient_analysis import TransientAnalysisResult
        

def solution(data_sources=None, channel=None):
    """Return a Result object which can provide information on a given set, on a given scoping...
    
    Parameters
    ----------
    Can be a filepath to the file you want to open, or a dpf.core.DataSources().
        
    Examples
    --------
    result = post.result("file.rst")
    """
    _model = Model(data_sources, channel)
    analysis_type = _model.metadata.result_info.analysis_type
    data_sources = _model.metadata.data_sources
    if(analysis_type == _AnalysisType.static):
        return StaticAnalysisResult(data_sources, _model)
    elif (analysis_type == _AnalysisType.modal):
        return ModalAnalysisResult(data_sources, _model)
    elif (analysis_type == _AnalysisType.harmonic):
        return HarmonicAnalysisResult(data_sources, _model)
    elif (analysis_type == _AnalysisType.transient):
        return TransientAnalysisResult(data_sources, _model)
    else:
        raise Exception("Unknown model.metadata.result_info.")
     
        
def build_docs(path=None):
    """Build HTML documentation.  This outputs all available
    operator types for the loaded operators.

    HTML is saved in the directory mentionned as parameters.
    
    Example
    -----
    from ansys.dpf import post
    post.build_doc("d:/temp/dpf_doc.html")
    
    Parameters
    -----
    str: output path for the documentation. Default is current directory.
    """
    if path is None:
        path = os.getcwd()
        path += "/dataProcessingDoc.html"
    doc_op = dpf.core.Operator('html_doc')
    doc_op.inputs.output_path.connect(path)
    doc_op.run()