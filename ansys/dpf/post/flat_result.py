##########################################################################
#                                                                        #
#          Copyright (C) 2020 ANSYS Inc.  All Rights Reserved            #
#                                                                        #
# This file contains proprietary software licensed from ANSYS Inc.       #
# This header must remain in any source code despite modifications or    #
# enhancements by any party.                                             #
#                                                                        #
##########################################################################
# Version: 1.0                                                           #
# Author(s): C.Bellot/R.Lagha/L.Paradis                                  #
# contact(s): ramdane.lagha@ansys.com                                    #
##########################################################################

"""Module containing the DpfFlatResult class and its childs. 
Each class highlights an analysis type, and provides hardcoded 
methods to get a result object regarding the wanted result type."""

from ansys.dpf.post.result_data import ResultData
from ansys.dpf.core import Operator as _Operator
from ansys.dpf.core.common import types

def _get_result_data_function_of_operator(name, instance, data_sources, **kwargs):
    location = None
    element_scoping = None
    node_scoping = None
    named_selection = None
    el_shape = None
    time_step = None
    grouping = None
    phase = None
    if "phase" in kwargs:
        if not isinstance(instance, DpfFlatComplexResult):
            raise Exception("Phase key-word argument can be used when the analysis types implies complex result (Harmonic analysis, Modal analysis...).")
        phase = kwargs["phase"]
    if "location" in kwargs:
        location = kwargs["location"]
    if "element_scoping" in kwargs:
        element_scoping = kwargs["element_scoping"]
    if "node_scoping" in kwargs:
        node_scoping = kwargs["node_scoping"]
    if "named_selection" in kwargs:
        named_selection = kwargs["named_selection"]
    if "el_shape" in kwargs:
        el_shape = kwargs["el_shape"]
    if "time_step" in kwargs:
        time_step = kwargs["time_step"]
    if "grouping" in kwargs:
        grouping = kwargs["grouping"]
    return ResultData(name, data_sources, location=location, element_scoping=element_scoping, 
             node_scoping = node_scoping, named_selection = named_selection, el_shape = el_shape, 
             time_step = time_step, grouping = grouping, phase = phase)


def _check_elemental_location(**kwargs):
    if "location" in kwargs:
        if(kwargs["location"] != 'Elemental'):
            raise Exception("Only an Elemental location can be ued with an elemental result.")

def _check_nodal_location(**kwargs):
    if "location" in kwargs:
        if(kwargs["location"] != 'Nodal'):
            raise Exception("Only a Nodal location can be ued with a nodal result.")

class DpfFlatResult:
    """Main class of flat result.
    
    Parameters
    ----
    None
    """
    def __init__(self, data_sources, metadata):
        self._data_sources = data_sources
        self._model_metadata = metadata
           
    def get_result_info(self):
        print(self._model_metadata.result_info)
    
    #nodal results
    def nodal_displacement(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("U", self, self._data_sources, **kwargs)
    def nodal_speed(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("V", self, self._data_sources, **kwargs)
    def nodal_acceleration(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("A", self, self._data_sources, **kwargs)
    def nodal_reaction_force(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("RF", self, self._data_sources, **kwargs)
    def nodal_force(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("F", self, self._data_sources, **kwargs)
    def nodal_moment(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("M", self, self._data_sources, **kwargs)
    def nodal_temperature(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("TEMP", self, self._data_sources, **kwargs)
    def nodal_raw_displacement(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("UTOT", self, self._data_sources, **kwargs)
    def nodal_raw_reaction_force(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("RFTOT", self, self._data_sources, **kwargs)
    def nodal_electric_field(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("EF", self, self._data_sources, **kwargs)
    def nodal_electric_potential(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("VOLT", self, self._data_sources, **kwargs)
    
    def modal_basis(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ModalBasis", self, self._data_sources, **kwargs)
    
    #element nodal results (get result at nodes or at elements)
    def elemental_stress(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("S", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_stress(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("S", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_elastic_strain(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_elastic_strain(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("EPEL", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_plastic_strain(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_plastic_strain(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("EPPL", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_structural_temperature(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("BFE", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_structural_temperature(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("BFE", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_thermal_strains(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ETH", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_thermal_strains(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ETH", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_eqv_stress_parameter(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_eqv_stress_parameter(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_SEPL", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_stress_ratio(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_stress_ratio(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_SRAT", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_hydrostatic_pressure(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_hydrostatic_pressure(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_HPRES", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_accu_eqv_plastic_strain(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_accu_eqv_plastic_strain(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_EPEQ", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_plastic_state_variable(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_plastic_state_variable(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_PSV", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_accu_eqv_creep_strain(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_accu_eqv_creep_strain(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_CREQ", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_plastic_strain_energy_density(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_plastic_strain_energy_density(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_PLWK", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_creep_strain_energy_density(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_creep_strain_energy_density(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_CRWK", self, self._data_sources, location="Nodal", **kwargs)
    def elemental_elastic_strain_energy_density(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Elemental", **kwargs)
    def nodal_elastic_strain_energy_density(self, **kwargs):
        _check_nodal_location(**kwargs)
        return _get_result_data_function_of_operator("ENL_ELENG", self, self._data_sources, location="Nodal", **kwargs)
    
    
    #element nodal result (here only elemental result given because the nodal one is given in nodal category)
    def elemental_force(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENF", self, self._data_sources, location="Elemental", **kwargs)
    
    #element results
    def elemental_contact_status(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_STAT", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_penetration(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_PENE", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_pressure(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_PRES", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_friction_stress(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_SFRIC", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_total_stress(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_STOT", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_sliding_distance(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_SLIDE", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contactgap_distance(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_GAP", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_surface_heat_flux(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_FLUX", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_contact_fluid_penetration_pressure(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ECT_FRES", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_volume(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_VOL", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_mass(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ElementalMass", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_stiffness_matrix_energy(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_SE", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_artificial_hourglass_energy(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_AHO", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_kinetic_energy(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_KE", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_co_energy(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_CO", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_incremental_energy(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_INC", self, self._data_sources, location="Elemental", **kwargs)
    def elemental_thermal_dissipation_energy(self, **kwargs):
        _check_elemental_location(**kwargs)
        return _get_result_data_function_of_operator("ENG_TH", self, self._data_sources, location="Elemental", **kwargs)

    
class DpfFlatComplexResult(DpfFlatResult):
    """Main class of flat result if the analysis gives complex result (Modal, Harmonic).
    
    Parameters
    ----
    None
    """
    def _get_amplitude_evaluation(self, result_data):
        # resultData = _get_result_data_function_of_operator(name, self, self._data_sources, **kwargs)
        resultData = result_data
        modulus_op = _Operator("modulus")
        modulus_op.inputs.fields_container.connect(resultData.result_fields_container)     
        resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        return resultData
    
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
    
class StaticAnalysisResult(DpfFlatResult):
    """Mechanic result class, which will provide all the API for Mechanic analysis."""
    
class ModalAnalysisResult(DpfFlatComplexResult):
    """Modal result class, which will provide all the API for Modal analysis."""
    
class HarmonicAnalysisResult(DpfFlatComplexResult):
    """Thermal result class, which will provide all the API for Thermal analysis with complex result."""
    
class TransientAnalysisResult(DpfFlatResult):
    """Transient result class, which will provide all the API for Transient analysis."""