"""This module contains all miscellaneous results."""

from ansys.dpf.core import Operator as _Operator
from ansys.dpf.core import locations

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.post.result_data import ResultData


class Misc:
    """Contains miscellaneous results.

    In this class, the phase keyword argument is available while calling results.
    """

    def __init__(self, model, data_sources):
        """Initialize this class."""
        self._model = model
        self._data_sources = data_sources

    # tools
    def _get_result_data_function_of_operator(
        self, name, instance, data_sources, b_elem_average: bool = False, **kwargs
    ):
        """Get a ``ResultData`` object with all available keywords."""
        location = None
        element_scoping = None
        node_scoping = None
        named_selection = None
        time = None
        grouping = None
        phase = None
        subresult = None
        set = None
        mapdl_grouping = None
        time_scoping = None
        if _AvailableKeywords._phase in kwargs:
            # if not isinstance(instance, dpf_solution.DpfComplexSolution):
            #     raise Exception(
            #         "Phase keyword argument can be used when the "
            #         "analysis type implies a complex result "
            #         "(harmonic analysis, modal analysis...).")
            #     )
            phase = kwargs[_AvailableKeywords._phase]
        if _AvailableKeywords.location in kwargs:
            location = kwargs[_AvailableKeywords.location]
        if _AvailableKeywords.element_scoping in kwargs:
            element_scoping = kwargs[_AvailableKeywords.element_scoping]
        if _AvailableKeywords.node_scoping in kwargs:
            node_scoping = kwargs[_AvailableKeywords.node_scoping]
        if _AvailableKeywords.named_selection in kwargs:
            named_selection = kwargs[_AvailableKeywords.named_selection]
        if _AvailableKeywords.time in kwargs:
            time = kwargs[_AvailableKeywords.time]
        if _AvailableKeywords.set in kwargs:
            set = kwargs[_AvailableKeywords.set]
        if _AvailableKeywords.grouping in kwargs:
            grouping = kwargs[_AvailableKeywords.grouping]
        if _AvailableKeywords.mapdl_grouping in kwargs:
            mapdl_grouping = kwargs[_AvailableKeywords.mapdl_grouping]
        if _AvailableKeywords._subresult in kwargs:
            subresult = kwargs[_AvailableKeywords._subresult]
        if _AvailableKeywords.time_scoping in kwargs:
            time_scoping = kwargs[_AvailableKeywords.time_scoping]
        return ResultData(
            name,
            data_sources,
            self._model,
            b_elem_average,
            location=location,
            element_scoping=element_scoping,
            node_scoping=node_scoping,
            named_selection=named_selection,
            time=time,
            grouping=grouping,
            phase=phase,
            subresult=subresult,
            mapdl_grouping=mapdl_grouping,
            set=set,
            time_scoping=time_scoping,
        )

    def _check_elemental_location(self, **kwargs):
        """Check if the location keyword with an elemental value is set.

        If the location keyword is not set, an exception is raised.
        """
        if _AvailableKeywords.location in kwargs:
            if kwargs[_AvailableKeywords.location] != locations.elemental:
                raise Exception(
                    "Only an elemental location can be used with an elemental result."
                )

    def _check_elemnodal_location(self, **kwargs):
        """Check if the location keyword with an elemental value is set.

        If the location keyword is not set, an exception is raised.
        """
        if _AvailableKeywords.location in kwargs:
            if kwargs[_AvailableKeywords.location] != locations.elemntal_nodal:
                raise Exception(
                    "Only an elemental nodal location can be used with an elemental nodal result."
                )

    def _check_nodal_location(self, **kwargs):
        """Check if the location keyword with an elemental value is set.

        If the location keyword is not set, n exception is ratised.

        """
        if _AvailableKeywords.location in kwargs:
            if kwargs[_AvailableKeywords.location] != locations.nodal:
                raise Exception(
                    "Only a nodal location can be used with a nodal result."
                )


class MecanicMisc(Misc):
    """Provides the miscellaneous mecanic result."""

    # nodal results
    def nodal_displacement(self, **kwargs):
        """Get result data for the nodal displacement."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "U", self, self._data_sources, **kwargs
        )

    def nodal_velocity(self, **kwargs):
        """Get result data for the nodal velocity."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "V", self, self._data_sources, **kwargs
        )

    def nodal_acceleration(self, **kwargs):
        """Get result data for the nodal acceleration."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "A", self, self._data_sources, **kwargs
        )

    def nodal_reaction_force(self, **kwargs):
        """Get result data for the nodal reaction force."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "RF", self, self._data_sources, **kwargs
        )

    def nodal_force(self, **kwargs):
        """Get result data for the nodal force."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "F", self, self._data_sources, **kwargs
        )

    def nodal_moment(self, **kwargs):
        """Get result data for the nodal moment."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "M", self, self._data_sources, **kwargs
        )

    def nodal_raw_displacement(self, **kwargs):
        """Get result data for the nodal raw displacement."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "UTOT", self, self._data_sources, **kwargs
        )

    def nodal_raw_reaction_force(self, **kwargs):
        """Get result data for the nodal raw reaction force."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "RFTOT", self, self._data_sources, **kwargs
        )

    def modal_basis(self, **kwargs):
        """Get result data for the modal basis."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ModalBasis", self, self._data_sources, **kwargs
        )

    # element nodal results (get result at nodes or at elements)
    def elemental_stress(self, **kwargs):
        """Get result data for the elemental stress."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "S",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_stress(self, **kwargs):
        """Get result data for the elemental nodal stress."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "S", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_stress(self, **kwargs):
        """Get result data for the nodal stress."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "S", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_elastic_strain(self, **kwargs):
        """Get result data for the elemental elastic strain."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EPEL",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_elastic_strain(self, **kwargs):
        """Get result data for the elemental nodal elastic strain."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EPEL", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_elastic_strain(self, **kwargs):
        """Get result data for the nodal elastic strain."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EPEL", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_plastic_strain(self, **kwargs):
        """Get result data for the elemental plastic strain."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EPPL",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_plastic_strain(self, **kwargs):
        """Get result data for the elemental nodal plastic strain."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EPPL", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_plastic_strain(self, **kwargs):
        """Get result data for the nodal plastic strain."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EPPL", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_structural_temperature(self, **kwargs):
        """Get result data for the elemental structural temperature."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "BFE",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_structural_temperature(self, **kwargs):
        """Get result data for the elemental nodal structural temperature."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "BFE", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_structural_temperature(self, **kwargs):
        """Get result data for the nodal structural temperature."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "BFE", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_thermal_strains(self, **kwargs):
        """Get result data for the elemental thermal strains."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ETH",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_thermal_strains(self, **kwargs):
        """Get result data for the elemental nodal thermal strains."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ETH", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_thermal_strains(self, **kwargs):
        """Get result data for the nodal thermal strains."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ETH", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_eqv_stress_parameter(self, **kwargs):
        """Get result data for an elemental equivalent stress parameter."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_SEPL",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_eqv_stress_parameter(self, **kwargs):
        """Get result data for the elemental nodal equivalent stress parameter."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_SEPL", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_eqv_stress_parameter(self, **kwargs):
        """Get result data for the nodal equivalent stress parameter."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_SEPL", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_stress_ratio(self, **kwargs):
        """Get result data for the elemental stress ratio."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_SRAT",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_stress_ratio(self, **kwargs):
        """Get result data for the elemental nodal stress ratio."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_SRAT", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_stress_ratio(self, **kwargs):
        """Get result data for the nodal stress ratio."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_SRAT", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_hydrostatic_pressure(self, **kwargs):
        """Get result data for the elemental hydrostatic pressure."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_HPRES",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_hydrostatic_pressure(self, **kwargs):
        """Get result data for the elemental nodal hydrostatic pressure."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_HPRES", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_hydrostatic_pressure(self, **kwargs):
        """Get result data for the nodal hydrostatic pressure."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_HPRES", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_accu_eqv_plastic_strain(self, **kwargs):
        """Get result data for the elemental accurate eqvivalent plastic strain."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_EPEQ",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_accu_eqv_plastic_strain(self, **kwargs):
        """Get result data for the elemental nodal accurate eqvivalent plastic strain."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_EPEQ", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_accu_eqv_plastic_strain(self, **kwargs):
        """Get result data for the nodal accurate equivalent plastic strain result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_EPEQ", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_plastic_state_variable(self, **kwargs):
        """Get result data for the elemental plastic state variable."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_PSV",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_plastic_state_variable(self, **kwargs):
        """Get result data for the elemental nodal plastic state variable."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_PSV", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_plastic_state_variable(self, **kwargs):
        """Get result data for the nodal plastic state variable."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_PSV", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_accu_eqv_creep_strain(self, **kwargs):
        """Get result data for the elemental accurate equivalent creep strain."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_CREQ",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_accu_eqv_creep_strain(self, **kwargs):
        """Get result data for the elemental nodal accurate equivalent creep strain."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_CREQ", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_accu_eqv_creep_strain(self, **kwargs):
        """Get result data for the nodal accurate equivalent creep strain."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_CREQ", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_plastic_strain_energy_density(self, **kwargs):
        """Get result data for the elemental plastic strain energy density."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_PLWK",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_plastic_strain_energy_density(self, **kwargs):
        """Get result data for the elemental nodal plastic strain energy density."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_PLWK", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_plastic_strain_energy_density(self, **kwargs):
        """Get result data for the nodal plastic strain energy density."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_PLWK", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_creep_strain_energy_density(self, **kwargs):
        """Get result data for the elemental creep strain energy density."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_CRWK",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_creep_strain_energy_density(self, **kwargs):
        """Get result data for the elemental nodal creep strain energy density."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_CRWK", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_creep_strain_energy_density(self, **kwargs):
        """Get result data for the nodal creep strain energy density."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_CRWK", self, self._data_sources, location="Nodal", **kwargs
        )

    def elemental_elastic_strain_energy_density(self, **kwargs):
        """Get result data for the elemental elastic strain energy density."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_ELENG",
            self,
            self._data_sources,
            location="Elemental",
            b_elem_average=True,
            **kwargs
        )

    def elemental_nodal_elastic_strain_energy_density(self, **kwargs):
        """Get result data for the elemental nodal elastic strain energy density."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_ELENG", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    def nodal_elastic_strain_energy_density(self, **kwargs):
        """Get result data for the nodal elastic strain energy density."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENL_ELENG", self, self._data_sources, location="Nodal", **kwargs
        )

    # element nodal result (here only elemental result given because
    # the nodal one is given in nodal category)
    def elemental_force(self, **kwargs):
        """Get result data for the elemental force."""
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator(
            "ENF", self, self._data_sources, location="Elemental", **kwargs
        )
        return self._elemental_nodal_to_elemental_result(resData)

    def elemental_nodal_force(self, **kwargs):
        """Get result data for the elemental nodal force."""
        self._check_elemnodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENF", self, self._data_sources, location="ElementalNodal", **kwargs
        )

    # element results
    def elemental_contact_status(self, **kwargs):
        """Get result data for the elemental contact status."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_STAT", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_penetration(self, **kwargs):
        """Get result data for the elemental contact penetration."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_PENE", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_pressure(self, **kwargs):
        """Get result data for the elemental contact pressure."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_PRES", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_friction_stress(self, **kwargs):
        """Get result data for the elemental contact friction stress."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_SFRIC", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_total_stress(self, **kwargs):
        """Get result data for the elemental contact total stress."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_STOT", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_sliding_distance(self, **kwargs):
        """Get result data for the elemental contact sliding distance."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_SLIDE", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contactgap_distance(self, **kwargs):
        """Get result data for the elemental contact gap distance."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_GAP", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_surface_heat_flux(self, **kwargs):
        """Get result data for the elemental contact surface heat flux."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_FLUX", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_contact_fluid_penetration_pressure(self, **kwargs):
        """Get result data for the elemental contact fluid penetration pressure."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ECT_FRES", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_volume(self, **kwargs):
        """Get result data for the elemental volume."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_VOL", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_mass(self, **kwargs):
        """Get result data for the elemental mass."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ElementalMass", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_stiffness_matrix_energy(self, **kwargs):
        """Get result data for the elemental stiffness matrix energy."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_SE", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_artificial_hourglass_energy(self, **kwargs):
        """Get result data for the elemental artificial hourglass energy."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_AHO", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_kinetic_energy(self, **kwargs):
        """Get result data for the elemental kinetic energy."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_KE", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_co_energy(self, **kwargs):
        """Get result data for the elemental co energy."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_CO", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_incremental_energy(self, **kwargs):
        """Get result data for the elemental incremental energy."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_INC", self, self._data_sources, location="Elemental", **kwargs
        )

    def elemental_thermal_dissipation_energy(self, **kwargs):
        """Get result data for the elemental thermal dissipation energy."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "ENG_TH", self, self._data_sources, location="Elemental", **kwargs
        )

    # special results
    def von_mises_stress(self, **kwargs):
        """Get result data for the nodal von Mises stress.

        The default location of this result is nodal.  Use the location keyword
        ``"ElementalNodal"`` to get an elemental nodal result or ``"Elemental"``
        to get an elemental result.
        """
        return self._get_result_data_function_of_operator(
            "S_eqv", self, self._data_sources, **kwargs
        )


class ComplexMecanicMisc(MecanicMisc):
    """Contains miscellaneous results."""

    # tools
    def _get_amplitude_evaluation(self, result_data):
        # resultData = self._get_result_data_function_of_operator(
        #     name, self, self._data_sources, **kwargs
        # )
        resultData = result_data
        modulus_op = _Operator("modulus")
        modulus_op.inputs.fields_container.connect(
            resultData._evaluator._result_operator.outputs.fields_container
        )
        resultData._evaluator._chained_operators[modulus_op.name] = (
            """This operator computes the amplitude of the result """
            """(when result has complex values)."""
        )
        # resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        resultData._evaluator._result_operator = modulus_op
        return resultData

    # results
    #!TODO
    def nodal_displacement_amplitude(self, **kwargs):
        """Get result data for the nodal displacement amplitude."""
        result_data = self.nodal_displacement(**kwargs)
        return self._get_amplitude_evaluation(result_data)

    def elemental_stress_amplitude(self, **kwargs):
        """Get result data for the elemental stress amplitude."""
        result_data = self.elemental_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)

    def nodal_stress_amplitude(self, **kwargs):
        """Get result data for the nodal stress amplitude."""
        result_data = self.nodal_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)


class ThermalMisc(Misc):
    """Contains miscellaneous results for thermal analysis.

    In this class the phase keyword argument is also available while calling results.
    """

    def nodal_electric_field(self, **kwargs):
        """Get result data for the nodal electric field."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "EF", self, self._data_sources, **kwargs
        )

    def nodal_electric_potential(self, **kwargs):
        """Get result data for the nodal electric potential."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "VOLT", self, self._data_sources, **kwargs
        )

    def nodal_temperature(self, **kwargs):
        """Get result data for the nodal thermal."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator(
            "TEMP", self, self._data_sources, **kwargs
        )
