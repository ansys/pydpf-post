"""Module containing the ``DpfSolution`` class and its children.

Each class highlights an analysis type and provides hard-coded methods to get
a result object regarding the wanted result type.

This module also contains the ``DpfComplexSolution`` class, which is a
child of the ``DpfSolution`` class.  In addition to the classic APIs, the complex
result introduces an amplitude evaluation.
"""

import re

from ansys.dpf.core import locations

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.post.displacement import ComplexDisplacement, Displacement
from ansys.dpf.post.electric_results import ElectricField, ElectricPotential
from ansys.dpf.post.errors import NodalLocationError
from ansys.dpf.post.misc_results import ComplexMecanicMisc, MecanicMisc
from ansys.dpf.post.strain import (
    ComplexElasticStrain,
    ComplexPlasticStrain,
    ElasticStrain,
    PlasticStrain,
)
from ansys.dpf.post.stress import ComplexStress, Stress
from ansys.dpf.post.temperature import (
    ComplexStructuralTemperature,
    HeatFlux,
    StructuralTemperature,
    Temperature,
)


class DpfSolution:
    """Provides the main class of the DPF-Post solution."""

    def __init__(self, data_sources, model):
        """Initialize the solution using ``data_sources`` and ``dpf.core.Model`` objects."""
        self._data_sources = data_sources
        self._model = model

    @property
    def mesh(self):
        """Mesh representation of the model.

        Returns the :class:`ansys.dpf.core.MeshedRegion` class.
        """
        return self._model.metadata.meshed_region

    @property
    def time_freq_support(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._model.metadata.time_freq_support

    def get_result_info(self):
        """Get result file information.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> print(solution.get_result_info()) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        Static analysis
        Unit system: MKS: m, kg, N, s, V, A, degC
        Physics Type: ...
        Available results:
             -  displacement: Nodal Displacement
             -  reaction_force: Nodal Force
             -  stress: ElementalNodal Stress
             -  elemental_volume: Elemental Volume
             -  stiffness_matrix_energy: Elemental Energy-stiffness matrix
             -  artificial_hourglass_energy: Elemental Hourglass Energy
             -  thermal_dissipation_energy: Elemental thermal dissipation energy
             -  kinetic_energy: Elemental Kinetic Energy
             -  co_energy: Elemental co-energy
             -  incremental_energy: Elemental incremental energy
             -  elastic_strain: ElementalNodal Strain
             -  structural_temperature: ElementalNodal Temperature
        """
        return self._model.metadata.result_info

    @staticmethod
    def _check_nodal_location(**kwargs):
        if _AvailableKeywords.location in kwargs:
            location = kwargs.pop(_AvailableKeywords.location)
            if location != locations.nodal:
                raise NodalLocationError()

    def __str__(self):
        """Get the string representation of this class."""
        txt = (
            "%s object." % re.sub(r"(?<!^)(?=[A-Z])", " ", type(self).__name__)
            + "\n\n\nData Sources\n------------------------------\n"
        )
        ds_str = self._data_sources.__str__()
        txt += ds_str
        txt += "\n\n"
        txt += self._model.__str__()
        return txt


class DpfMecanicSolution(DpfSolution):
    """Provides the mecanic solution."""

    def __init__(self, data_sources, model):
        """Initialize this class."""
        super().__init__(data_sources, model)
        self.misc = MecanicMisc(model, data_sources)

    # result classes
    def stress(self, **kwargs):
        """Get a stress object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress(node_scoping = [1, 43])
        """
        return Stress(data_sources=self._data_sources, model=self._model, **kwargs)

    def elastic_strain(self, **kwargs):
        """Get an elastic strain object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> elastic_strain = solution.elastic_strain(node_scoping = [1, 43])
        """
        return ElasticStrain(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def plastic_strain(self, **kwargs):
        """Return a plastic strain object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> plastic_strain = solution.plastic_strain(node_scoping = [1, 43])
        """
        return PlasticStrain(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def displacement(self, **kwargs):
        """Get a displacement object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> displacement = solution.displacement(node_scoping = [1, 43])
        """
        return Displacement(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def structural_temperature(self, **kwargs):
        """Get a temperature object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
        >>> temperature = solution.structural_temperature(node_scoping = [1, 43])
        """
        return StructuralTemperature(
            data_sources=self._data_sources, model=self._model, **kwargs
        )


class DpfMecanicComplexSolution(DpfSolution):
    """Provides the main class of the DPF-Post solution if the analysis gives a complex solution."""

    def __init__(self, data_sources, model):
        """Initialize this class."""
        super().__init__(data_sources, model)
        self.misc = ComplexMecanicMisc(model, data_sources)

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nThis may contain complex results."

    def has_complex_result(self):
        """Check if the solution object has complex values (complex frequencies).

        Returns
        -------
        bool
            ``True`` if the solution has complex results, ``False`` otherwise.
        """
        tfq_sup = self._model.metadata.time_freq_support
        if not tfq_sup:
            return False
        if tfq_sup.complex_frequencies is None:
            return False
        return True

    def displacement(self, **kwargs):
        """Get a displacement object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> displacement = solution.displacement(node_scoping = [1, 43])
        """
        return ComplexDisplacement(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def structural_temperature(self, **kwargs):
        """Get a temperature object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
        >>> temperature = solution.structural_temperature(node_scoping = [1, 43])
        """
        return ComplexStructuralTemperature(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def plastic_strain(self, **kwargs):
        """Get a plastic strain object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> plastic_strain = solution.plastic_strain(node_scoping = [1, 43])
        """
        return ComplexPlasticStrain(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def elastic_strain(self, **kwargs):
        """Get an elastic strain object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> elastic_strain = solution.elastic_strain(node_scoping = [1, 43])
        """
        return ComplexElasticStrain(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def stress(self, **kwargs):
        """Get a stress object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress(node_scoping = [1, 43])
        """
        return ComplexStress(
            data_sources=self._data_sources, model=self._model, **kwargs
        )


class DpfThermalSolution(DpfSolution):
    """Provides the main class of the DPF-Post solution for a thermal analysis type."""

    def __init__(self, data_sources, model):
        """Initialize this class."""
        super().__init__(data_sources, model)
        # self.misc = ThermalMisc(model, data_sources)

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nThis may contain complex results."

    def temperature(self, **kwargs):
        """Get a temperature object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.steady_therm)
        >>> temp = solution.temperature(node_scoping = [1, 43])
        """
        self._check_nodal_location(**kwargs)
        return Temperature(data_sources=self._data_sources, model=self._model, **kwargs)

    def heat_flux(self, **kwargs):
        """Get a heat flux object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.steady_therm)
        >>> hf = solution.heat_flux(node_scoping = [1, 43])
        """
        return HeatFlux(data_sources=self._data_sources, model=self._model, **kwargs)

    def electric_potential(self, **kwargs):
        """Get an electric potential object, from which you can possibly get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.electric_therm)
        >>> ep = solution.electric_potential(node_scoping = [1, 43])
        """
        self._check_nodal_location(**kwargs)
        return ElectricPotential(
            data_sources=self._data_sources, model=self._model, **kwargs
        )

    def electric_field(self, **kwargs):
        """Get an electric field object from which you can possible to get result data.

        Parameters
        ----------
        **kwargs
            List of keyword arguments. You can use the :class:`print_available_keywords
            <ansys.dpf.post.print_available_keywords>` method to find keyword arguments.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.electric_therm)
        >>> ef = solution.electric_field(node_scoping = [1, 43])
        """
        return ElectricField(
            data_sources=self._data_sources, model=self._model, **kwargs
        )
