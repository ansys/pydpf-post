"""Module containing the method to instantiate the result object.

Post-utility
------------

This module is used for the initialization of PyDPF-Post objects.
"""
from typing import TypeVar, Union
import warnings

from ansys.dpf.core.model import Model
from ansys.dpf.core.server_types import BaseServer

from ansys.dpf.post.common import (
    AvailableSimulationTypes,
    _AnalysisType,
    _AvailableKeywords,
    _PhysicsType,
    simulation_type_str_to_class,
)
from ansys.dpf.post.harmonic_analysis import HarmonicAnalysisSolution
from ansys.dpf.post.modal_analysis import ModalAnalysisSolution
from ansys.dpf.post.simulation import Simulation
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

    .. deprecated:: 3.0
       Use :func:`load_simulation` instead.

    This class provides information on a given set on a given scoping.

    Parameters
    ----------
    data_sources: str, ansys.dpf.core.DataSources
         Path to the file to open or the :class:`ansys.dpf.core.DataSources` class.
    physics_type: common._PhysicsType, str, optional
        Type of phsyics described in the specified data sources. Options are
        ``"mecanic"`` or ``"thermal"``. The default is ``None``, in which case
        the data sources are read to determine the physics type.
    analysis_type: common._AnalysisType, str, optional
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
                Warning(
                    "Physics type is defaulting to 'mechanical'. Specify the physics type",
                    "keyword if it is invalid.",
                )
            )
            physics_type = _PhysicsType.mechanical

    if not analysis_type:
        try:
            analysis_type = _model.metadata.result_info.analysis_type
        except Exception as e:
            warnings.warn(
                Warning(
                    "Analysis type is defaulting to 'static'. Specify the analysis"
                    "type keyword if it is invalid."
                )
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


SimulationType = TypeVar("SimulationType", bound=Simulation)


def load_simulation(
    data_sources,
    simulation_type: Union[AvailableSimulationTypes, str] = None,
    server: Union[BaseServer, None] = None,
) -> SimulationType:
    """Loads a simulation and returns a :class:`ansys.dpf.post.simulation.Simulation` object.

    This class provides the main interface to explore and manipulate results, meshes, geometries,
    and other entities associated with the result files given in input.
    The interface exposed depends on the type of simulation selected with the argument
    `simulation_type`. Each one proposes a post-processing context with its specific vocabulary and
    most common post-processing functionalities.
    Available simulation types are listed in
    :class:`<AvailableSimulationTypes> ansys.dpf.post.common.AvailableSimulationTypes`.

    Parameters
    ----------
    data_sources: str, ansys.dpf.core.DataSources
         Path to the file to open or the :class:`ansys.dpf.core.DataSources` class.
    simulation_type:
        Type of simulation to create when loading the specified data sources.
        Each type of simulation gives access to specific properties and methods to better fit
        expectations and vocabulary of each context.
        This defaults to the simulation type corresponding to the combination of physics type
        and analysis type detected automatically by DPF when reading the result files.
        If nothing is detected, this will default to a static mechanical type of simulation.
        The best practice is to define this parameter to select the right post-processing context.
        Options are given in
        :class:`<AvailableSimulationTypes> ansys.dpf.post.common.AvailableSimulationTypes`.
    server:
        DPF server connected to a remote or local instance.

    Returns
    -------
    An instance of one of the subclasses of the
    :class:`Simulation <ansys.dpf.post.simulation.Simulation>` class.

    .. versionadded:: 3.0
        This function replaces the deprecated :func:`load_solution` function.

    Examples
    --------
    Load the example static result.

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> simulation = post.load_simulation(examples.static_rst)
    """
    _model = Model(data_sources, server=server)
    data_sources = _model.metadata.data_sources

    if not simulation_type:
        try:
            physics_type = _model.metadata.result_info.physics_type
        except Exception as e:
            warnings.warn(
                Warning(
                    "Physics type is defaulting to 'mechanical'. Specify 'simulation_type' "
                    "keyword if you want to use another type."
                )
            )
            physics_type = _PhysicsType.mechanical
        try:
            analysis_type = _model.metadata.result_info.analysis_type
        except Exception as e:
            warnings.warn(
                Warning(
                    "Analysis type is set to 'static' as default. Specify the 'simulation_type' "
                    "keyword if you want to use another type."
                )
            )
            analysis_type = _AnalysisType.static

        if physics_type == _PhysicsType.thermal:
            if analysis_type == _AnalysisType.static:
                raise NotImplementedError
            elif analysis_type == _AnalysisType.transient:
                raise NotImplementedError
            else:
                raise ValueError(
                    f"Unknown analysis type '{analysis_type}' for {physics_type}."
                )
        elif (
            physics_type == _PhysicsType.mecanic
            or physics_type == _PhysicsType.mechanical
        ):
            if analysis_type == _AnalysisType.static:
                simulation_type = AvailableSimulationTypes.static_mechanical
            elif analysis_type == _AnalysisType.modal:
                simulation_type = AvailableSimulationTypes.modal_mechanical
            elif (
                analysis_type == _AnalysisType.harmonic
                or analysis_type == _AnalysisType.msup
            ):
                simulation_type = AvailableSimulationTypes.harmonic_mechanical
            elif analysis_type == _AnalysisType.transient:
                simulation_type = AvailableSimulationTypes.transient_mechanical
            else:
                raise ValueError(
                    f"Unknown analysis type '{analysis_type}' for {physics_type}."
                )
        elif physics_type == _PhysicsType.fluid:
            if analysis_type in [_AnalysisType.steady, _AnalysisType.static]:
                simulation_type = AvailableSimulationTypes.steady_fluid
            elif analysis_type in [_AnalysisType.unsteady, _AnalysisType.transient]:
                simulation_type = AvailableSimulationTypes.unsteady_fluid
            else:
                raise ValueError(
                    f"Unknown analysis type '{analysis_type}' for {physics_type}."
                )
        else:
            raise ValueError(f"Unknown physics type '{physics_type}.")

    if simulation_type in [
        getattr(AvailableSimulationTypes, x) for x in vars(AvailableSimulationTypes)
    ]:
        return simulation_type_str_to_class[simulation_type](
            data_sources, server=_model._server
        )
    else:
        raise ValueError(
            f"Simulation type '{simulation_type}' is not a recognized simulation type."
            f" Please use ansys.dpf.post.common.AvailableSimulationTypes or instantiate one of the"
            f" available Simulation sub-classes directly."
        )


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
