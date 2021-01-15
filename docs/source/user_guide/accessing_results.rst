.. _user_guide_accessing_results:

*****************
Accessing results
*****************

The DPF-Post module provides an simplified Python interface to DPF, 
thus enabling rapid post-processing. To proceed, a solution object 
must be instantiated first.
Refer to the previous part to proceed (:ref:`_user_guide_post_processing'). 

Once the solution is instantiated, it is possible to access to the 
contained results.

Results can me manipulated as result object's instances. The following 
code shows how to proceed:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate a result object 
	>>> # here is called a displacement result
	>>> displacement = solution.displacement()
	>>> # stress, elastic_strain (...) can also be called. 
	>>> # Refer to the list below to know which result object can be instantiated.
	
**Each result object can be instantiated with a bunch of keyword arguments.** 
For a **full example using keyword arguments**, see :ref:`_ref_result_keywords`.
	
A result object is available for each main DPF-Post result type. 
See the following list to know which result can be accessed this way. 

List of available result objects using DPF-Post API:
- displacement
- stress
- elastic_strain
- plastic_strain
- structural_temperature

**Before calling a result object, you need to be sure that the result information 
is contained is your result file.** 
	

Displacement
------------
Displacement can be called as a result object, using the following code:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the displacement result object 
	>>> displacement = solution.displacement()

Its **location** can **only be nodal**. 

Then, the data can be accessed through a result data object. The whole ResultData 
API can be found following :ref:`_ref_result_data`.

For example, if the subresult "UY" is wanted (displacement Y component), the 
following proceedure can be followed:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the displacement result object 
	>>> displacement = solution.displacement()
	>>> # get the y displacement result data 
	>>> u_y = displacement.y
	>>> u_y.get_data_at_field()

Other subresults and the whole vector data can be accessed the same way.

.. autoclass:: ansys.dpf.post.displacement.Displacement
    :members:
	

Stress
------
Stress can be called as a result object, using the following code:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the stress result object 
	>>> stress = solution.stress()
	
Then, the data can be accessed through a result data object. The whole ResultData 
API can be found following :ref:`_ref_result_data`.

For example, if the subresult "SYY" is wanted (stress YY component), the 
following proceedure can be followed:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the displacement result object 
	>>> stress = solution.stress()
	>>> # get the yy stress result data 
	>>> s_yy = stress.yy
	>>> s_yy.get_data_at_field()

Other subresults and the whole tensor data can be accessed the same way.

.. autoclass:: ansys.dpf.post.stress.Stress
    :members:
	

Strain (elastic, plastic)
-------------------------
Elastic or plastic strain can be called as a result object, using the following code:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the strain result objects
	>>> elastic_strain = solution.elastic_strain()
	>>> plastic_strain = solution.plastic_strain()
	
Then, the data can be accessed through a result data object. The whole ResultData 
API can be found following :ref:`_ref_result_data`.

For example, if the elastic strain XY component is wanted, the 
following proceedure can be followed:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the displacement result object 
	>>> elastic_strain = solution.elastic_strain()
	>>> # get the xy elastic strain result data 
	>>> e_yy = elastic_strain.xy
	>>> e_yy.get_data_at_field()

Other subresults and the whole tensor data can be accessed the same way.

.. autoclass:: ansys.dpf.post.strain.ElasticStrain
    :members:
	
.. autoclass:: ansys.dpf.post.strain.PlasticStrain
    :members:	

Structural temperature
----------------------
Structural temperature (or system temperature) can be called as a result object, using the following code:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the elastic strain result object 
	>>> structural_temperature = solution.structural_temperature()
	
Then, the data can be accessed through a result data object. The whole ResultData 
API can be found following :ref:`_ref_result_data`.

To access the temperature scalar data, the following proceedure can be followed:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the displacement result object 
	>>> structural_temperature = solution.structural_temperature()
	>>> # get the xy elastic strain result data 
	>>> temperature = structural_temperature.scalar
	>>> temperature.get_data_at_field()

.. autoclass:: ansys.dpf.post.temperature.StructuralTemperature
    :members:
