"""Module containing the specific Transient Analysis Solution class."""

from ansys.dpf.post.dpf_solution import DpfMechanicSolution, DpfThermalSolution


class TransientAnalysisSolution(DpfMechanicSolution):
    """Transient solution class, which will provide all the API for Transient analysis."""


class ThermalTransientAnalysisSolution(DpfThermalSolution):
    """Transient solution class, which will provide all the API for Transient analysis."""
