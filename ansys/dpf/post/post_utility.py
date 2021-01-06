"""Module containing the method to instantiate the result
object. Initialization of post objects.
"""

from ansys.dpf.core.model import Model
 
from ansys.dpf.post.common import _AnalysisType, _AvailableKeywords, _PhysicsType
from ansys.dpf.post.static_analysis import (StaticAnalysisSolution,
                                            ThermalStaticAnalysisSolution)
from ansys.dpf.post.modal_analysis import ModalAnalysisSolution
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.transient_analysis import (TransientAnalysisSolution,
                                               ThermalTransientAnalysisSolution)


def load_solution(data_sources):
    """Return a ``Result`` object which can provide information on a given
    set, on a given scoping.

    Parameters
    ----------
    data_sources : str or dpf.core.DataSources
         filepath to the file you want to open, or a dpf.core.DataSources().

    Examples
    --------
    >>> solution = post.solution("file.rst")
    """
    _model = Model(data_sources)
    data_sources = _model.metadata.data_sources

    analysis_type = _model.metadata.result_info.analysis_type
    physics_type = _model.metadata.result_info.physics_type
    if physics_type == _PhysicsType.thermal:
        if analysis_type == _AnalysisType.static:
            return ThermalStaticAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.transient:
            return ThermalTransientAnalysisSolution(data_sources, _model)
        else:
            raise Exception(f"Unknown analysis type '{analysis_type}' for thermal.")
    elif physics_type == _PhysicsType.mecanic:
        if analysis_type == _AnalysisType.static:
            return StaticAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.modal:
            return ModalAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.harmonic:
            return HarmonicAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.transient:
            return TransientAnalysisSolution(data_sources, _model)
        else:
            raise Exception(f"Unknown analysis type '{analysis_type}' for mechanical.")
    else:
        raise Exception(f"Unknown physics type '{physics_type}.")


def print_available_keywords():
    """Print the keywords that can be used into the result object.

    Examples
    --------
    >>> from ansys.dpf import post
    >>> solution = post.load_solution('file.rst')
    >>> stress = solution.stress(PRINTED_KEYWORDS_CAN_be_USED_HERE)
    """
    txt = _AvailableKeywords().__str__()
    print(txt)
