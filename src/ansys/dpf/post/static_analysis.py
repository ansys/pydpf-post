"""Module containing the specific Static Analysis Solution class."""

from ansys.dpf.post.dpf_solution import DpfMechanicSolution, DpfThermalSolution


class StaticAnalysisSolution(DpfMechanicSolution):
    """Static solution class, which will provide all the API for Static analysis."""


class ThermalStaticAnalysisSolution(DpfThermalSolution):
    """Static solution class, which will provide all the API for Static analysis."""
