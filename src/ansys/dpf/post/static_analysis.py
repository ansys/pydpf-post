"""Module containing the specific Static Analysis Solution class."""

from ansys.dpf.post.dpf_solution import DpfMecanicSolution, DpfThermalSolution


class StaticAnalysisSolution(DpfMecanicSolution):
    """Static solution class, which will provide all the API for Static analysis."""


class ThermalStaticAnalysisSolution(DpfThermalSolution):
    """Static solution class, which will provide all the API for Static analysis."""
