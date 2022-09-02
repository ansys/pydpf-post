"""Module containing the method to instantiate the result object.

This module is used for the initialization of DPF-Post objects.
"""
from ansys.dpf.core.model import Model

from ansys.dpf.post.common import _AnalysisType, _AvailableKeywords, _PhysicsType
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.modal_analysis import ModalAnalysisSolution
from ansys.dpf.post.static_analysis import (
    StaticAnalysisSolution,
    ThermalStaticAnalysisSolution,
)
from ansys.dpf.post.transient_analysis import (
    ThermalTransientAnalysisSolution,
    TransientAnalysisSolution,
)


def load_solution(data_sources, physics_type=None, analysis_type=None):
    """Loads a solution and returns a :class:`ansys.dpf.post.Result` object.

    This class provides information on a given set on a given scoping.

    Parameters
    ----------
    data_sources : str or ansys.dpf.core.DataSources
         Path to the file to open or the :class:`ansys.dpf.core.DataSources` class.
    physics_type : common._PhysicsType, str, optional
        Type of phsyics described in the specified data sources. Options are
        ``"mecanic"`` or ``"thermal"``. The default is ``None``, in which case
        the data sources are read to determine the physics type.
    analysis_type : common._AnalysisType, str, optional
        Type of analysis described in the specified data sources. Options are
        ``"static"``, ``"modal"``, ``"harmonic"``, and ``"transient"``. The
        default is ``None``, in which case the data sources are read to determine
        the analysis type.

    Examples
    --------
    Load the example static result.

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.static_rst)
    """
    _model = Model(data_sources)
    data_sources = _model.metadata.data_sources

    if not physics_type:
        try:
            physics_type = _model.metadata.result_info.physics_type
        except Exception as e:
            warnings.warn(
                "Physics type is defaulting to 'mecanic'. Specify the physics type",
                "keyword if it is invalid.",
            )
            physics_type = _PhysicsType.mecanic

    if not analysis_type:
        try:
            analysis_type = _model.metadata.result_info.analysis_type
        except Exception as e:
            warnings.warn(
                "Analysis type is defaulting to 'static'. Specify the analysis"
                "type keyword if it is invalid.",
            )
            analysis_type = _AnalysisType.static

    if physics_type == _PhysicsType.thermal:
        if analysis_type == _AnalysisType.static:
            return ThermalStaticAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.transient:
            return ThermalTransientAnalysisSolution(data_sources, _model)
        else:
            raise ValueError(f"Unknown analysis type '{analysis_type}' for thermal.")
    elif (
        physics_type == _PhysicsType.mecanic or physics_type == _PhysicsType.mechanical
    ):
        if analysis_type == _AnalysisType.static:
            return StaticAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.modal:
            return ModalAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.harmonic:
            return HarmonicAnalysisSolution(data_sources, _model)
        elif analysis_type == _AnalysisType.transient:
            return TransientAnalysisSolution(data_sources, _model)
        else:
            raise ValueError(f"Unknown analysis type '{analysis_type}' for mechanical.")
    else:
        raise ValueError(f"Unknown physics type '{physics_type}.")


def print_available_keywords():
    """Print the keywords that can be used into the result object.

    Examples
    --------
    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
    >>> post.print_available_keywords() # doctest: +NORMALIZE_WHITESPACE
    element_scoping: list, int or dpf.core.Scoping
    grouping: str. Use post.grouping.(...) as helper.
    location: str. Use post.locations.(...) as helper.
    mapdl_grouping: int. Write 186 to get mapdl_elements solid_186.
    named_selection: str. Name of named_selection.
    node_scoping: list, int or dpf.core.Scoping
    path: DpfPath object that
                contains a list of coordinates,
                e.g. [[0.1, 0.0, 0.0],
                [0.0, 0.1, 0.0]].
    set: int
    time: float
    time_scoping: list, int or dpf.core.Scoping
    """
    txt = _AvailableKeywords().__str__()
    print(txt)
