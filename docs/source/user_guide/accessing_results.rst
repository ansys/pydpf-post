.. _user_guide_accessing_results:

*****************
Accessing results
*****************

The ``Solution`` object is the entry point for browsing the contents of a 
result file (see :ref:`user_guide_post_processing`). It also provides access 
to the results themselves. These are contained in ``Result`` objects that can be
returned from dedicated methods:  

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
	
A number of options, formulated as **keyword arguments**, can be used to further specify the result 
type, scope, and time, among others. For detailed examples see :ref:`ref_result_keywords`.

At present, DPF-Post supports two types of result files: 

* Structural (\*.rst) 
* Thermal/Electric (\*.rth)
	
Once loaded into the ``Solution`` instance, a result file will offer its own variety of ``Result`` 
objects in accordance with its type.

**Only request results types that are available in the file.** To verify their presence
see :ref:`user_guide_accessing_file_metadata`.

Structural (\*.rst) result files
================================

Upon loading a ``Solution`` object from a structural analysis results file (\*.rst), the following
``Result`` objects may be queried: 

* displacement
* stress
* elastic_strain
* plastic_strain
* structural_temperature

Displacement
------------
Displacement is the DOF solution for structural analyses.
Its ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the displacement result object 
	>>> displacement = solution.displacement()

Being a DOF solution, its location argument **must be nodal** (inferred above). 

The displacement ``Result`` corresponds to a vector field. To obtain the scalar components of this field, 
e.g., y-components, access the following subresult:

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

For more details, see :ref:`ref_api_result_data`.

.. autoclass:: ansys.dpf.post.displacement.Displacement
	:inherited-members:
	

Stress
------
The stress ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the stress result object 
	>>> stress = solution.stress()
	
The stress ``Result`` corresponds to a tensor field. To obtain the scalar components of this field, 
e.g., the normal y-stresses, access the following subresult: 

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

Other components, as well as whole tensor data, can be queried accordingly.
For more details, see :ref:`ref_api_result_data`.

.. autoclass:: ansys.dpf.post.stress.Stress
	:inherited-members:
	

Strain (elastic, plastic)
-------------------------
The elastic and plastic strain ``Result`` objects can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the strain result objects
	>>> elastic_strain = solution.elastic_strain()
	>>> plastic_strain = solution.plastic_strain()
	
A strain ``Result`` corresponds to a tensor field. To obtain the scalar components of this field, 
e.g., shear xy-strains, access the following subresult:

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

Other components, as well as whole tensor data, can be queried accordingly.
For more details, see :ref:`ref_api_result_data`.

.. autoclass:: ansys.dpf.post.strain.ElasticStrain
	:inherited-members:
	
.. autoclass:: ansys.dpf.post.strain.PlasticStrain
	:inherited-members:

Structural temperature
----------------------
The structural (or system) temperature ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the elastic strain result object 
	>>> structural_temperature = solution.structural_temperature()
	
To access the temperature scalar field use the following:

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
	:inherited-members:
	
	
Miscellaneous results
---------------------

Other ``Result`` objects may be available in the ``Solution``: 

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # get the result data
	>>> acceleration = solution.misc.nodal_acceleration()

All **keyword arguments** are available for miscellaneous results (see :ref:`ref_result_keywords`),
except location. Notably, some subresults may be available as keyword arguments,
for example, the scalar components of acceleration:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # get the result data
	>>> acceleration = solution.misc.nodal_acceleration(subresult="Y")  
	
Verify the result file contents to determine available queries (see :ref:`user_guide_accessing_file_metadata`).

.. autoclass:: ansys.dpf.post.misc_results.MecanicMisc
	:members:

Thermal/Electric (\*.rth) result files
======================================

Upon loading a ``Solution`` object from a thermal/electric analysis results file (\*.rth), the following
``Result`` objects may be queried: 

* temperature
* heat flux
* electric field
* electric potential

Temperature
-----------
Temperature is the DOF solution for thermal analyses.  Its ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.steady_therm)
	>>> # instantiate the temperature result object 
	>>> temperature = solution.temperature()

Being a DOF solution, its location argument **must be nodal** (inferred above).  The scalar field can 
be obtained directly:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the temperature result object 
	>>> temperature = solution.temperature()
	>>> # get the y temperature result data 
	>>> temp = temperature.scalar
	>>> temp.get_data_at_field()

.. autoclass:: ansys.dpf.post.temperature.Temperature
	:inherited-members:


Heat flux
---------
The heat flux ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.steady_therm)
	>>> # instantiate the heat_flux result object 
	>>> heat_flux = solution.heat_flux()

The heat flux ``Result`` corresponds to a vector field. To obtain the scalar components of this field, 
e.g., x-components, access the following subresult:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the temperature result object 
	>>> heat_flux = solution.heat_flux()
	>>> # get the y heat_flux result data 
	>>> heat_flux_x = heat_flux.x
	>>> heat_flux_x.get_data_at_field()
	
Other components can be queried accordingly.
For more details, see :ref:`ref_api_result_data`.

.. autoclass:: ansys.dpf.post.temperature.HeatFlux
	:inherited-members:
	

Electric field
--------------
The electric field ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.electric_therm)
	>>> # instantiate the electric field result object 
	>>> electric_field = solution.electric_field()

The electric field ``Result`` corresponds to a vector field.  To obtain the scalar components of this field,
e.g., x-components, acces the following subresult:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the electric field result object 
	>>> electric_field = solution.electric_field()
	>>> # get the y electricfield result data 
	>>> electric_field_x = electric_field.x
	>>> electric_field_x.get_data_at_field()
	
For more details, see :ref:`ref_api_result_data`.

.. autoclass:: ansys.dpf.post.electric_results.ElectricField
	:inherited-members:
	
	
Electric potential
------------------
The electric potential ``Result`` object can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.steady_therm)
	>>> # instantiate the electric potential result object 
	>>> electric_potential = solution.electric_potential()

It represents a scalar field whose values can be obtained as follows:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate the electric potential result object 
	>>> electric_potential = solution.electric_potential()
	>>> # get the y electric potential result data 
	>>> ep = electric_potential.scalar
	>>> ep.get_data_at_field()

.. autoclass:: ansys.dpf.post.electric_results.ElectricPotential
	:inherited-members: