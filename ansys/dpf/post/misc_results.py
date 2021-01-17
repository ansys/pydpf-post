"""This module contains all the miscallenous results."""

from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.core import locations
from ansys.dpf.post import dpf_solution
from ansys.dpf.core import Operator as _Operator

class Misc():
    """This class contains miscellaneous results.
    
    Here the phase keyword is also available while calling results.
    """
    
    def __init__(self, model, data_sources):
        self._model = model
        self._data_sources = data_sources
    
    #tools
    def _get_result_data_function_of_operator(self, name, instance, data_sources, b_elem_average: bool = False, **kwargs):
        """This method check which are the used keywords, then retusn a ResultData instance 
        computed with all available keywords."""
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
            #     raise Exception("Phase key-word argument can be used when the analysis types implies complex result (Harmonic analysis, Modal analysis...).")
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
        return ResultData(name, data_sources, self._model, b_elem_average, location=location, element_scoping=element_scoping, 
                 node_scoping = node_scoping, named_selection = named_selection,
                 time = time, grouping = grouping, phase = phase, subresult=subresult, mapdl_grouping=mapdl_grouping, 
                 set=set, time_scoping=time_scoping)
    
    def _check_elemental_location(self, **kwargs):
        """Check if the location keyword with an Elemental value is set. If not, raise Exception."""
        if _AvailableKeywords.location in kwargs:
            if(kwargs[_AvailableKeywords.location] != locations.elemental):
                raise Exception("Only an Elemental location can be used with an elemental result.")

    def _check_nodal_location(self, **kwargs):
        """Check if the location keyword with an Elemental value is set. If not, raise Exception."""
        if _AvailableKeywords.location in kwargs:
            if(kwargs[_AvailableKeywords.location] != locations.nodal):
                raise Exception("Only a Nodal location can be used with a nodal result.")
    
    
class MecanicMisc(Misc):
    
    #nodal results
    def nodal_displacement(self, **kwargs):
        """Returns a nodal_displacement result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("U", self, self._data_sources, **kwargs)
    
    def nodal_velocity(self, **kwargs):
        """Returns a nodal_velocity result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("V", self, self._data_sources, **kwargs)
    
    def nodal_acceleration(self, **kwargs):
        """Returns a nodal_acceleration result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("A", self, self._data_sources, **kwargs)
    
    def nodal_reaction_force(self, **kwargs):
        """Returns a nodal_reaction_force result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("RF", self, self._data_sources, **kwargs)
    
    def nodal_force(self, **kwargs):
        """Returns a nodal_force result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("F", self, self._data_sources, **kwargs)
    
    def nodal_moment(self, **kwargs):
        """Returns a nodal_moment result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("M", self, self._data_sources, **kwargs)
    
    def nodal_raw_displacement(self, **kwargs):
        """Returns a nodal_raw_displacement result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("UTOT", self, self._data_sources, **kwargs)
    
    def nodal_raw_reaction_force(self, **kwargs):
        """Returns a nodal_raw_reaction_force result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("RFTOT", self, self._data_sources, **kwargs)
    
    def modal_basis(self, **kwargs):
        """Returns a modal_basis result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ModalBasis", self, self._data_sources, **kwargs)
    
    
    #element nodal results (get result at nodes or at elements)
    def elemental_stress(self, **kwargs):
        """Returns a elemental_stress result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("S", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_stress(self, **kwargs):
        """Returns a elemental_nodal_stress result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("S", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_stress(self, **kwargs):
        """Returns a nodal_stress result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("S", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_elastic_strain(self, **kwargs):
        """Returns a elemental_elastic_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_elastic_strain(self, **kwargs):
        """Returns a elemental_nodal_elastic_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_elastic_strain(self, **kwargs):
        """Returns a nodal_elastic_strain result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_plastic_strain(self, **kwargs):
        """Returns a elemental_plastic_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_plastic_strain(self, **kwargs):
        """Returns a elemental_nodal_plastic_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_plastic_strain(self, **kwargs):
        """Returns a nodal_plastic_strain result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_structural_temperature(self, **kwargs):
        """Returns a elemental_structural_temperature result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("BFE", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_structural_temperature(self, **kwargs):
        """Returns a elemental_nodal_structural_temperature result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("BFE", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_structural_temperature(self, **kwargs):
        """Returns a nodal_structural_temperature result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("BFE", self, self._data_sources, location="Nodal", **kwargs)
        
    def elemental_thermal_strains(self, **kwargs):
        """Returns a elemental_thermal_strains result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ETH", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_thermal_strains(self, **kwargs):
        """Returns a elemental_nodal_thermal_strains result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ETH", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_thermal_strains(self, **kwargs):
        """Returns a nodal_thermal_strains result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ETH", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_eqv_stress_parameter(self, **kwargs):
        """Returns a elemental_eqv_stress_parameter result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_eqv_stress_parameter(self, **kwargs):
        """Returns a elemental_nodal_eqv_stress_parameter result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_eqv_stress_parameter(self, **kwargs):
        """Returns a nodal_eqv_stress_parameter result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_stress_ratio(self, **kwargs):
        """Returns a elemental_stress_ratio result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_stress_ratio(self, **kwargs):
        """Returns a elemental_nodal_stress_ratio result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_stress_ratio(self, **kwargs):
        """Returns a nodal_stress_ratio result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_hydrostatic_pressure(self, **kwargs):
        """Returns a elemental_hydrostatic_pressure result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_hydrostatic_pressure(self, **kwargs):
        """Returns a elemental_nodal_hydrostatic_pressure result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_hydrostatic_pressure(self, **kwargs):
        """Returns a nodal_hydrostatic_pressure result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_accu_eqv_plastic_strain(self, **kwargs):
        """Returns a elemental_accu_eqv_plastic_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_accu_eqv_plastic_strain(self, **kwargs):
        """Returns a elemental_nodal_accu_eqv_plastic_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_accu_eqv_plastic_strain(self, **kwargs):
        """Returns a nodal_accu_eqv_plastic_strain result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_plastic_state_variable(self, **kwargs):
        """Returns a elemental_plastic_state_variable result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_plastic_state_variable(self, **kwargs):
        """Returns a elemental_nodal_plastic_state_variable result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_plastic_state_variable(self, **kwargs):
        """Returns a nodal_plastic_state_variable result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_accu_eqv_creep_strain(self, **kwargs):
        """Returns a elemental_accu_eqv_creep_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_accu_eqv_creep_strain(self, **kwargs):
        """Returns a elemental_nodal_accu_eqv_creep_strain result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_accu_eqv_creep_strain(self, **kwargs):
        """Returns a nodal_accu_eqv_creep_strain result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_plastic_strain_energy_density(self, **kwargs):
        """Returns a elemental_plastic_strain_energy_density result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_plastic_strain_energy_density(self, **kwargs):
        """Returns a elemental_nodal_plastic_strain_energy_density result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_plastic_strain_energy_density(self, **kwargs):
        """Returns a nodal_plastic_strain_energy_density result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_creep_strain_energy_density(self, **kwargs):
        """Returns a elemental_creep_strain_energy_density result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_creep_strain_energy_density(self, **kwargs):
        """Returns a elemental_nodal_creep_strain_energy_density result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_creep_strain_energy_density(self, **kwargs):
        """Returns a nodal_creep_strain_energy_density result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_elastic_strain_energy_density(self, **kwargs):
        """Returns a elemental_elastic_strain_energy_density result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Elemental", b_elem_average = True, **kwargs)
    
    def elemental_nodal_elastic_strain_energy_density(self, **kwargs):
        """Returns a elemental_nodal_elastic_strain_energy_density result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_elastic_strain_energy_density(self, **kwargs):
        """Returns a nodal_elastic_strain_energy_density result data."""
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Nodal", **kwargs)
    
    
    #element nodal result (here only elemental result given because the nodal one is given in nodal category)
    def elemental_force(self, **kwargs):
        """Returns a elemental_force result data."""
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENF", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_force(self, **kwargs):
        """Returns a elemental_nodal_force result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENF", self, self._data_sources, location="Elemental", **kwargs)
    
    
    #element results
    def elemental_contact_status(self, **kwargs):
        """Returns a elemental_contact_status result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_STAT", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_penetration(self, **kwargs):
        """Returns a elemental_contact_penetration result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_PENE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_pressure(self, **kwargs):
        """Returns a elemental_contact_pressure result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_PRES", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_friction_stress(self, **kwargs):
        """Returns a elemental_contact_friction_stress result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_SFRIC", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_total_stress(self, **kwargs):
        """Returns a elemental_contact_total_stress result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_STOT", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_sliding_distance(self, **kwargs):
        """Returns a elemental_contact_sliding_distance result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_SLIDE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contactgap_distance(self, **kwargs):
        """Returns a elemental_contactgap_distance result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_GAP", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_surface_heat_flux(self, **kwargs):
        """Returns a elemental_contact_surface_heat_flux result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_FLUX", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_fluid_penetration_pressure(self, **kwargs):
        """Returns a elemental_contact_fluid_penetration_pressure result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_FRES", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_volume(self, **kwargs):
        """Returns a elemental_volume result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_VOL", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_mass(self, **kwargs):
        """Returns a elemental_mass result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ElementalMass", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_stiffness_matrix_energy(self, **kwargs):
        """Returns a elemental_stiffness_matrix_energy result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_SE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_artificial_hourglass_energy(self, **kwargs):
        """Returns a elemental_artificial_hourglass_energy result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_AHO", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_kinetic_energy(self, **kwargs):
        """Returns a elemental_kinetic_energy result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_KE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_co_energy(self, **kwargs):
        """Returns a elemental_co_energy result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_CO", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_incremental_energy(self, **kwargs):
        """Returns a elemental_incremental_energy result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_INC", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_thermal_dissipation_energy(self, **kwargs):
        """Returns a elemental_thermal_dissipation_energy result data."""
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_TH", self, self._data_sources, location="Elemental", **kwargs)


    #special results
    def von_mises_stress(self, **kwargs):
        """von_mises_stress output default location is Nodal. 
        Use the location keyword with 'ElementalNodal' value to 
        get an ElementalNodal result, and the 'Elemental' value 
        to get Elemental result.
        """
        return self._get_result_data_function_of_operator("S_eqv", self, self._data_sources, **kwargs)


class ComplexMecanicMisc(MecanicMisc):
    """This class contains miscellaneous results."""
    
    #tools
    def _get_amplitude_evaluation(self, result_data):
        # resultData = self._get_result_data_function_of_operator(name, self, self._data_sources, **kwargs)
        resultData = result_data
        modulus_op = _Operator("modulus")
        modulus_op.inputs.fields_container.connect(resultData._evaluator._result_operator.outputs.fields_container) 
        resultData._evaluator._chained_operators[modulus_op.name] = """This operator will compute the amplitude of the result (when result has complex values)."""
        # resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        resultData._evaluator._result_operator = modulus_op
        return resultData
    
    #results
    #!TODO
    def nodal_displacement_amplitude(self, **kwargs):
        result_data = self.nodal_displacement(**kwargs)
        return self._get_amplitude_evaluation(result_data)
    
    def elemental_stress_amplitude(self, **kwargs):
        result_data = self.elemental_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)
    
    def nodal_stress_amplitude(self, **kwargs):
        result_data = self.nodal_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)  
    
    
class ThermalMisc(Misc):
    """This class contains miscellaneous results for thermal analysis.
    
    Here the phase keyword is also available while calling results.
    """
    
    def nodal_electric_field(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("EF", self, self._data_sources, **kwargs)
    
    def nodal_electric_potential(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("VOLT", self, self._data_sources, **kwargs)
    
    def nodal_temperature(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("TEMP", self, self._data_sources, **kwargs)
    
    