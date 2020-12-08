"""Module containing the DpfSolution class and its childs. 
Each class highlights an analysis type, and provides hardcoded 
methods to get a result object regarding the wanted result type.

Module containing also the DpfComplexSolution class, child of DpfSolution class.
Additionnaly to the classic APIs, the complex result introduces an amplitude evaluation."""


from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.core import Operator as _Operator

from ansys.dpf.post.stress import Stress, ComplexStress
from ansys.dpf.post.strain import ElasticStrain, ComplexElasticStrain
from ansys.dpf.post.strain import PlasticStrain, ComplexPlasticStrain
from ansys.dpf.post.temperature import Temperature, ComplexTemperature
from ansys.dpf.post.displacement import Displacement, ComplexDisplacement

from ansys import dpf
        

class DpfSolution:
    """Main class of post result API.
    
    Parameters
    ----------
    None
    """
    def __init__(self, data_sources, model):
        self._data_sources = data_sources
        self._model = model
        self.misc = dpf.post.misc_results.Misc(self._data_sources, self._model)
           
    def get_result_info(self):
        """Returns information about the result file.
        
        Parameters
        ----------
        None
        
        Examples
        --------
        The following code:
        >>> from ansys.dpf import post
        >>> result = post.result("file.rst")
        >>> print(result.get_result_info())
        
        Will return:
            Static analysis
            Unit system: Metric (m, kg, N, s, V, A)
            Physics Type: Mecanic
            Available results:
                 -  displacement                                        
                 -  volume                                        
        """
        return self._model.metadata.result_info
    
    def __str__(self):
        txt = '%s result object.' % self._model.metadata.result_info.analysis_type.capitalize() +\
        '\n\n\nData Sources\n------------------------------\n'
        ds_str = self._data_sources.__str__()
        txt += ds_str
        txt += "\n\n"
        txt += self._model.__str__()
        return txt
        
    
    #tools
    def _check_phase(self, **kwargs):
        if _AvailableKeywords.phase in kwargs:
            if not isinstance(self, DpfComplexSolution):
                raise Exception("Phase key-word argument can be used when the analysis types implies complex result (Harmonic analysis, Modal analysis...).")
    
    #result classes            
    def stress(self, **kwargs):
        self._check_phase(**kwargs)
        return Stress(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def elastic_strain(self, **kwargs):
        self._check_phase(**kwargs)
        return ElasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def plastic_strain(self, **kwargs):
        self._check_phase(**kwargs)
        return PlasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def displacement(self, **kwargs):
        self._check_phase(**kwargs)
        return Displacement(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def temperature(self, **kwargs):
        self._check_phase(**kwargs)
        return Temperature(data_sources=self._data_sources, model=self._model, **kwargs)
    

class DpfComplexSolution(DpfSolution):
    """Main class of post solution if the analysis gives complex solution (Modal, Harmonic).
    
    Parameters
    ----------
    None
    """
    def _get_amplitude_evaluation(self, result_data):
        # resultData = self._get_result_data_function_of_operator(name, self, self._data_sources, **kwargs)
        resultData = result_data
        modulus_op = _Operator("modulus")
        modulus_op.inputs.fields_container.connect(resultData._evaluator._result_operator.outputs.fields_container) 
        resultData._evaluator._chained_operators[modulus_op.name] = """This operator will compute the amplitude of the result (when result has complex values)."""
        # resultData.result_fields_container = modulus_op.get_output(0, types.fields_container)
        resultData._evaluator._result_operator = modulus_op
        return resultData
    
    def __str__(self):
        txt = super().__str__()
        txt += "\n"
        txt += "This can contain complex result."
        return txt
    
    def has_complex_result(self):
        """Tests if the solution object has complex values (check the complex frequencies).
        
        Returns
        -------
        Boolean (True if has_complex_result)"""
        tfq_sup = self._model.metadata.time_freq_support
        if not tfq_sup:
            return False
        if (tfq_sup.complex_frequencies == None):
            return False
        return True
    
    
    def complex_displacement(self, **kwargs):
        self._check_phase(**kwargs)
        return ComplexDisplacement(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def complex_temperature(self, **kwargs):
        self._check_phase(**kwargs)
        return ComplexTemperature(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def complex_plastic_strain(self, **kwargs):
        self._check_phase(**kwargs)
        return ComplexPlasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def complex_elastic_strain(self, **kwargs):
        self._check_phase(**kwargs)
        return ComplexElasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def complex_stress(self, **kwargs):
        self._check_phase(**kwargs)
        return ComplexStress(data_sources=self._data_sources, model=self._model, **kwargs)
    
    
    def _nodal_displacement_amplitude(self, **kwargs):
        result_data = self._nodal_displacement(**kwargs)
        return self._get_amplitude_evaluation(result_data)
    
    def _elemental_stress_amplitude(self, **kwargs):
        result_data = self._elemental_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)
    
    def _nodal_stress_amplitude(self, **kwargs):
        result_data = self._nodal_stress(**kwargs)
        return self._get_amplitude_evaluation(result_data)    

    