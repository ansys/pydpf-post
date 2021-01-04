.. _dpf_model_functions_ref:

API Reference
=============
Details of the DPF API.


Solution object
---------------
The solution object instantiates an object that is built on the result file. 
Use the following code to instantiate a solution object.

.. code:: python

    >>> from ansys.dpf import post
	>>> solution = post.load_solution(r'../../../tests/testfiles/model_with_ns.rst')


Result object
-------------
The result object can be manipulated to get different result data. 

.. code:: python

    >>> from ansys.dpf import post
	>>> solution = post.load_solution(r'../../../tests/testfiles/model_with_ns.rst')
	>>> # Displacement result object
	>>> displacement = solution.displacement()
	>>> # Stress result object
	>>> stress = solution.stress()
	>>> # Elastic strain result object
	>>> elastic_strain = solution.elastic_strain()


ResultData Class
---------------------
The ResultData class is instantiated from a result object. 
It enables an easy access to the data. The following code shows 
how to get a ResultData instance. 
	
.. code:: python

    >>> from ansys.dpf import post
	>>> solution = post.load_solution(r'../../../testfiles/tests/model_with_ns.rst')
	>>> displacement = solution.displacement()
	>>> result = displacement.vector
	>>> # access the data
	>>> vector.get_data_at_field(0)
	
.. autoclass:: ansys.dpf.post.result_data.ResultData
    :members:
	


