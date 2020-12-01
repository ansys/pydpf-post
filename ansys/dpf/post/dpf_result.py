"""Module containing the DpfResult class and its childs. 
Each class highlights an analysis type, and provides hardcoded 
methods to get a result object regarding the wanted result type.

Module containing also the DpfComplexResult class, child of DpfResult class.
Additionnaly to the classic APIs, the complex result introduces an amplitude evaluation."""


from ansys.dpf.post.result_data import ResultData
from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.core import Operator as _Operator
from ansys.dpf.core.common import locations
        

class DpfResult:
    """Main class of post result API.
    
    Parameters
    ----
    None
    """
    def __init__(self, data_sources, model):
        self._data_sources = data_sources
        self._model = model
           
    def get_result_info(self):
        """Returns information about the result file.
        
        Parameters
        -----
        None
        
        Example
        -----
        The following code:
            from ansys.dpf import post
            result = post.result("file.rst")
            print(result.get_result_info())
        
        Will return:
            Static analysis
            Unit system: Metric (m, kg, N, s, V, A)
            Physics Type: Mecanic
            Available results:
                 -  displacement                                        
                 -  volume                                        
        -----
        """
        return self._model.metadata.result_info
    
    def __str__(self):
        txt = '%s result object.\n\n' % self._model_metadata.result_info.analysis_type.capitalize() +\
        'The open dataSource has contains the following information:\n'
        ds_str = self._data_sources.__str__()
        txt += ds_str
        return txt
        
        
    #tools
    def _get_result_data_function_of_operator(self, name, instance, data_sources, **kwargs):
        """This method check which are the used keywords, then retusn a ResultData instance 
        computed with all available keywords."""
        location = None
        element_scoping = None
        node_scoping = None
        named_selection = None
        el_shape = None
        time = None
        grouping = None
        phase = None
        subresult = None
        set = None
        mapdl_grouping = None
        time_scoping = None
        if _AvailableKeywords.phase in kwargs:
            if not isinstance(instance, DpfComplexResult):
                raise Exception("Phase key-word argument can be used when the analysis types implies complex result (Harmonic analysis, Modal analysis...).")
            phase = kwargs[_AvailableKeywords.phase]
        if _AvailableKeywords.location in kwargs:
            location = kwargs[_AvailableKeywords.location]
        if _AvailableKeywords.element_scoping in kwargs:
            element_scoping = kwargs[_AvailableKeywords.element_scoping]
        if _AvailableKeywords.node_scoping in kwargs:
            node_scoping = kwargs[_AvailableKeywords.node_scoping]
        if _AvailableKeywords.named_selection in kwargs:
            named_selection = kwargs[_AvailableKeywords.named_selection]
        if _AvailableKeywords.el_shape in kwargs:
            el_shape = kwargs[_AvailableKeywords.el_shape]
        if _AvailableKeywords.time in kwargs:
            time = kwargs[_AvailableKeywords.time]
        if _AvailableKeywords.set in kwargs:
            set = kwargs[_AvailableKeywords.set]
        if _AvailableKeywords.grouping in kwargs:
            grouping = kwargs[_AvailableKeywords.grouping]
        if _AvailableKeywords.mapdl_grouping in kwargs:
            mapdl_grouping = kwargs[_AvailableKeywords.mapdl_grouping]
        if _AvailableKeywords.subresult in kwargs:
            subresult = kwargs[_AvailableKeywords.subresult]
        if _AvailableKeywords.time_scoping in kwargs:
            time_scoping = kwargs[_AvailableKeywords.time_scoping]
        return ResultData(name, data_sources, self._model, self, location=location, element_scoping=element_scoping, 
                 node_scoping = node_scoping, named_selection = named_selection, el_shape = el_shape, 
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
                
                
    def _elemental_nodal_to_elemental_result(self, result_data):
        avg = _Operator("to_elemental_fc")
        avg.inputs.fields_container.connect(result_data._result_operator.outputs.fields_container)
        result_data.result_fields_container = avg.outputs.fields_container()
        return result_data
        
    
    
    #nodal results
    def nodal_displacement(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("U", self, self._data_sources, **kwargs)
    
    def nodal_speed(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("V", self, self._data_sources, **kwargs)
    
    def nodal_acceleration(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("A", self, self._data_sources, **kwargs)
    
    def nodal_reaction_force(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("RF", self, self._data_sources, **kwargs)
    
    def nodal_force(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("F", self, self._data_sources, **kwargs)
    
    def nodal_moment(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("M", self, self._data_sources, **kwargs)
    
    def nodal_temperature(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("TEMP", self, self._data_sources, **kwargs)
    
    def nodal_raw_displacement(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("UTOT", self, self._data_sources, **kwargs)
    
    def nodal_raw_reaction_force(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("RFTOT", self, self._data_sources, **kwargs)
    
    def nodal_electric_field(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("EF", self, self._data_sources, **kwargs)
    
    def nodal_electric_potential(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("VOLT", self, self._data_sources, **kwargs)
    
    def modal_basis(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ModalBasis", self, self._data_sources, **kwargs)
    
    
    #element nodal results (get result at nodes or at elements)
    def elemental_stress(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("S", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_stress(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("S", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_stress(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("S", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_elastic_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_elastic_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_elastic_strain(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_plastic_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_plastic_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_plastic_strain(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_structural_temperature(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("BFE", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_structural_temperature(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("BFE", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_structural_temperature(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("BFE", self, self._data_sources, location="Nodal", **kwargs)
        
    def elemental_thermal_strains(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ETH", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_thermal_strains(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ETH", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_thermal_strains(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ETH", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_eqv_stress_parameter(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_eqv_stress_parameter(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_eqv_stress_parameter(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_stress_ratio(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_stress_ratio(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_stress_ratio(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_hydrostatic_pressure(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_hydrostatic_pressure(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_hydrostatic_pressure(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_accu_eqv_plastic_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_accu_eqv_plastic_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_accu_eqv_plastic_strain(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_plastic_state_variable(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_plastic_state_variable(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_plastic_state_variable(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_accu_eqv_creep_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_accu_eqv_creep_strain(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_accu_eqv_creep_strain(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_plastic_strain_energy_density(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_plastic_strain_energy_density(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_plastic_strain_energy_density(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_creep_strain_energy_density(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_creep_strain_energy_density(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_creep_strain_energy_density(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Nodal", **kwargs)
    
    def elemental_elastic_strain_energy_density(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_elastic_strain_energy_density(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Elemental", **kwargs)
    
    def nodal_elastic_strain_energy_density(self, **kwargs):
        self._check_nodal_location(**kwargs)
        return self._get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Nodal", **kwargs)
    
    
    #element nodal result (here only elemental result given because the nodal one is given in nodal category)
    def elemental_force(self, **kwargs):
        self._check_elemental_location(**kwargs)
        resData = self._get_result_data_function_of_operator("ENF", self, self._data_sources, location="Elemental", **kwargs)
        return self._elemental_nodal_to_elemental_result(resData)
    
    def elemental_nodal_force(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENF", self, self._data_sources, location="Elemental", **kwargs)
    
    
    #element results
    def elemental_contact_status(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_STAT", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_penetration(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_PENE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_pressure(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_PRES", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_friction_stress(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_SFRIC", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_total_stress(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_STOT", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_sliding_distance(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_SLIDE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contactgap_distance(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_GAP", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_surface_heat_flux(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_FLUX", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_contact_fluid_penetration_pressure(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ECT_FRES", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_volume(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_VOL", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_mass(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ElementalMass", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_stiffness_matrix_energy(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_SE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_artificial_hourglass_energy(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_AHO", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_kinetic_energy(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_KE", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_co_energy(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_CO", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_incremental_energy(self, **kwargs):
        self._check_elemental_location(**kwargs)
        return self._get_result_data_function_of_operator("ENG_INC", self, self._data_sources, location="Elemental", **kwargs)
    
    def elemental_thermal_dissipation_energy(self, **kwargs):
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


class DpfComplexResult(DpfResult):
    """Main class of post result if the analysis gives complex result (Modal, Harmonic).
    
    Parameters
    ----
    None
    """
    def _get_amplitude_evaluation(self, result_data):
        # resultData = self._get_result_data_function_of_operator(name, self, self._data_sources, **kwargs)
        resultData = result_data
        modulus_op = _Operator("modulus")
        modulus_op.inputs.fields_container.connect(resultData._result_operator.outputs.fields_container) 
        resultData._chained_operators[modulus_op.name] = """This operator will compute the amplitude of the result (when result has complex values)."""
        # resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        resultData._result_operator = modulus_op
        return resultData
    
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This can contain complex result."
        return txt
    
    def has_complex_result(self):
        """Tests if the result object has complex values (check the complex frequencies).
        
        Returns
        -----
        Boolean (True if has_complex_result)"""
        tfq_sup = self._model.metadata.time_freq_support
        if not tfq_sup:
            return False
        if (tfq_sup.complex_frequencies == None):
            return False
        return True
    
    
    #!TODO all results
    def nodal_displacement_amplitude(self, **kwargs):
        result_data = self.nodal_displacement(**kwargs)
        return self._get_amplitude_evaluation(result_data)
    
    def elemental_stress_amplitude(self, **kwargs):
        result_data = self.elemental_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)
    
    def nodal_stress_amplitude(self, **kwargs):
        result_data = self.nodal_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)    

    