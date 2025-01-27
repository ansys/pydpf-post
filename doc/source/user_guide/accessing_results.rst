.. _user_guide_accessing_results:

**************
Access results
**************

After a result file is loaded, you can access the results themselves. You can query results using dedicated methods.

This code shows how to instantiate the :class:`Simulation <ansys.dpf.post.simulation.Simulation>` object and then
get the ``displacement`` result for a loaded simulation:

.. code:: python

	Instantiate the simulation object

	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> simulation = post.load_simulation(examples.multishells_rst)

	Extract the displacement data as a DataFrame object

	>>> displacement = simulation.displacement()
	>>> # stress, elastic_strain (...) can also be called.

	See the following list for the results that can be
	extracted.

You can use *keyword arguments* to specify additional options,
including the components, scope, and time. For detailed examples,
see :ref:`ref_result_keywords`.

PyDPF-Post supports two types of result files:

* Structural (RST)
* Thermal/electric (RTH) (with the legacy ``load_solution()`` method only)

It also supports some LS-DYNA result files and some Fluent and CFX result files as indicated in
:ref:`solver_result_files_support`.



You should only request results available in the result file.
To determine which results are available, see :ref:`user_guide_accessing_file_metadata`.

Structural result files
=======================

You can use the legacy ``Solution`` object to access structural results.

After loading a ``Solution`` object from a structural analysis result (RST)
file, you can query these ``Result`` objects:

* ``displacement``
* ``stress``
* ``elastic_strain``
* ``plastic_strain``
* ``structural_temperature``

Displacement
------------
Displacement is the DOF solution for a structural analysis. The location argument
for a DOF solution must be modal.

This code shows how to access the :class:`Displacement <ansys.dpf.post.displacement.Displacement>`
result object:

.. code:: python

    # Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the displacement result object

    >>> displacement = solution.displacement()

The displacement ``Result`` object corresponds to a vector field. To obtain the scalar
components (y-components) of this field, access the subresult with this code:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the displacement result object

    >>> displacement = solution.displacement()

    Get the y displacement result data

    >>> u_y = displacement.y
    >>> u_y.get_data_at_field()

For more information, see :ref:`ref_api_result_data`.


Stress
------
This code shows to access the :class:`Stress <ansys.dpf.post.stress.Stress>` result
object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    # Instantiate the stress result object

    >>> stress = solution.stress()

A ``Stress`` result object corresponds to a tensor field. To obtain the scalar
components of this field, such as the normal y-stresses, access the subresult:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the stress result object

    >>> stress = solution.stress()

    Get the yy stress result data

    >>> s_yy = stress.yy
    >>> s_yy.get_data_at_field()

You can query other components, as well as whole tensor data, accordingly.
For more information, see :ref:`ref_api_result_data`.


Elastic and plastic strain
--------------------------
This code shows how to access the :class:`ElasticStrain <ansys.dpf.post.strain.ElasticStrain>` and
:class:`PlasticStrain <ansys.dpf.post.strain.PlasticStrain>` result objects:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the strain result objects

    >>> elastic_strain = solution.elastic_strain()
    >>> plastic_strain = solution.plastic_strain()

A ``Strain`` result object corresponds to a tensor field. To obtain the scalar
components of this field, such as the shear xy-strains, access the subresult:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the elastic strain result object

    >>> elastic_strain = solution.elastic_strain()

    Get the xy elastic strain result data

    >>> e_yy = elastic_strain.xy
    >>> e_yy.get_data_at_field()

You can query other components, as well as whole tensor data, accordingly.
For more information, see :ref:`ref_api_result_data`.


Structural temperature
----------------------
This code shows to access the :class:`StructuralTemperature <ansys.dpf.post.temperature.StructuralTemperature>`
result object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the structural temperature result object

    >>> structural_temperature = solution.structural_temperature()

To access the temperature scalar field, use this code:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the structural temperature result object

    >>> structural_temperature = solution.structural_temperature()

    Get the structural temperature result data

    >>> temperature = structural_temperature.scalar
    >>> temperature.get_data_at_field()


Miscellaneous results
---------------------
The ``Solution`` object might contain other miscellaneous :class:`ansys.dpf.post.misc_results.MecanicMisc`
result objects that you can access. For example, this code shows how to access the ``nodal_acceleration``
result object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Get the nodal acceleration result data

    >>> acceleration = solution.misc.nodal_acceleration()

All keyword arguments are available for miscellaneous results, except ``location``.
For more information, see :ref:`ref_result_keywords`.

Some subresults might be available as keyword arguments, such as the scalar
components of nodal acceleration:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Get the result data

    >>> acceleration = solution.misc.nodal_acceleration(subresult="Y")

To determine available queries, you can browse the metadata in the result file. For more
information, see :ref:`user_guide_accessing_file_metadata`.


Thermal/electric result files
=============================

You can access thermal/electric results using the legacy ``Solution`` object.

After loading a ``Solution`` object from a thermal/electric analysis
result file (RTH), you can query these ``Result`` objects:

* ``temperature``
* ``heat_flux``
* ``electric_field``
* ``electric_potential``

Temperature
-----------
Temperature is the DOF solution for a thermal analysis.

This code shows how to access the :class:`Temperature <ansys.dpf.post.temperature.Temperature>`
result object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.steady_therm)

    Instantiate the temperature result object

    >>> temperature = solution.temperature()

As inferred above, the ``location`` argument for a DOF solution must be nodal.
This code shows how to access the scalar field directly:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the temperature result object

    >>> temperature = solution.temperature()

    Get the y temperature result data

    >>> temp = temperature.scalar
    >>> temp.get_data_at_field()


Heat flux
---------
This code shows how to access the :class:`HeatFlux <ansys.dpf.post.temperature.HeatFlux>` result
object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.steady_therm)

    Instantiate the heat_flux result object

    >>> heat_flux = solution.heat_flux()


The ``HeatFlux`` result object corresponds to a vector field. To obtain the scalar
components (x-components) of this field, access the subresult:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the temperature result object

    >>> heat_flux = solution.heat_flux()

    Get the y heat_flux result data

    >>> heat_flux_x = heat_flux.x
    >>> heat_flux_x.get_data_at_field()

You can query other components accordingly. For more information, see
:ref:`ref_api_result_data`.


Electric field
--------------
This code shows how to access the :class:`ElectricField <ansys.dpf.post.electric_results.ElectricField>`
result object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.electric_therm)

    Instantiate the electric field result object

    >>> electric_field = solution.electric_field()

The ``electric_field`` result object corresponds to a vector field. To
obtain the scalar components of this field, such as the x-components, access
the subresult:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the electric field result object

    >>> electric_field = solution.electric_field()

    Get the y electricfield result data

    >>> electric_field_x = electric_field.x
    >>> electric_field_x.get_data_at_field()

For more information, see :ref:`ref_api_result_data`.


Electric potential
------------------
This code shows how to access the :class:`ElectricPotential <ansys.dpf.post.electric_results.ElectricPotential>`
result object:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.steady_therm)

    Instantiate the electric potential result object

    >>> electric_potential = solution.electric_potential()

The ``electric_potential`` result object corresponds to a scalar field. You can access
its values with:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the electric potential result object

    >>> electric_potential = solution.electric_potential()

    Get the y electric potential result data

    >>> ep = electric_potential.scalar
    >>> ep.get_data_at_field()
