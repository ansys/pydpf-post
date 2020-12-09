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

from ansys.dpf.post.misc_results import Misc, ComplexMisc
        

class DpfSolution:
    """Main class of post result API.
    
    Parameters
    ----------
    None
    """
    def __init__(self, data_sources, model):
        self._data_sources = data_sources
        self._model = model
        self.misc = Misc(model, data_sources)
           
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
        
    #result classes            
    def stress(self, **kwargs):
        """Returns a stress object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
            
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> stress = solution.stress(node_scoping = [1, 43])
        """
        return Stress(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def elastic_strain(self, **kwargs):
        """Returns an elastic strain object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> elastic_strain = solution.elastic_strain(node_scoping = [1, 43])
        """
        return ElasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def plastic_strain(self, **kwargs):
        """Returns a plastic strain object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> plastic_strain = solution.plastic_strain(node_scoping = [1, 43])
        """
        return PlasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def displacement(self, **kwargs):
        """Returns a displacement object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> displacement = solution.displacement(node_scoping = [1, 43])
        """
        return Displacement(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def temperature(self, **kwargs):
        """Returns a temperature object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> temperature = solution.temperature(node_scoping = [1, 43])
        """
        return Temperature(data_sources=self._data_sources, model=self._model, **kwargs)
    

class DpfComplexSolution(DpfSolution):
    """Main class of post solution if the analysis gives complex solution (Modal, Harmonic).
    
    Parameters
    ----------
    None
    """
    def __init__(self, data_sources, model):
        super().__init__(data_sources, model)
        self.misc = ComplexMisc(model, data_sources)
    
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
    
    
    def displacement(self, **kwargs):
        """Returns a displacement object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> displacement = solution.displacement(node_scoping = [1, 43])
        """
        return ComplexDisplacement(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def temperature(self, **kwargs):
        """Returns a temperature object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> temperature = solution.temperature(node_scoping = [1, 43])
        """
        return ComplexTemperature(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def plastic_strain(self, **kwargs):
        """Returns a plastic strain object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> plastic_strain = solution.plastic_strain(node_scoping = [1, 43])
        """
        return ComplexPlasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def elastic_strain(self, **kwargs):
        """Returns an elastic strain object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> elastic_strain = solution.elastic_strain(node_scoping = [1, 43])
        """
        return ComplexElasticStrain(data_sources=self._data_sources, model=self._model, **kwargs)
    
    def stress(self, **kwargs):
        """Returns a stress object from which it is possible to get ResultData.
        
        Parameters
        ----------
        **kwargs
            The list of keyword-arguments can be found using post.print_available_keywords().
        
        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.load_solution(file.rst)
        >>> stress = solution.stress(node_scoping = [1, 43])
        """
        return ComplexStress(data_sources=self._data_sources, model=self._model, **kwargs)
    
  

    