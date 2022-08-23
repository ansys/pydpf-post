.. _user_guide_accessing_results:

*****************
Accessing Results
*****************

The :class:`Solution <ansys.dpf.post.dpf_solution.DpfSolution>`
object is the entry point for browsing the contents of a result file
(see :ref:`user_guide_post_processing`). It also provides access to
the results themselves. These are contained in :class:`Result
<ansys.dpf.post.result_object.Result>` objects that can be returned from
dedicated methods:

.. code:: python

	Instantiate the solution object

	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)

	Instantiate a displacement result object

	>>> displacement = solution.displacement()
	>>> # stress, elastic_strain (...) can also be called.

	Refer to the list below to know which result object can be
	instantiated.

A number of options, formulated as **keyword arguments**, can be used
to further specify the result type, scope, and time, among others. For
detailed examples see :ref:`ref_result_keywords`.

At present, DPF-Post supports two types of result files:

* Structural ``(\*.rst)``
* Thermal/Electric ``(\*.rth)``

Once loaded into a :py:class:`Solution
<ansys.dpf.post.dpf_solution.DpfSolution>`, a result file will offer its own
variety of :class:`Result <ansys.dpf.post.result_object.Result>` objects in
accordance with its type.

**Only request results types that are available in the file.** To
verify their presence see :ref:`user_guide_accessing_file_metadata`.

Structural (\*.rst) result files
================================

Upon loading a ``Solution`` object from a structural analysis results
file (\*.rst), the following ``Result`` objects may be queried:

* ``displacement``
* ``stress``
* ``elastic_strain``
* ``plastic_strain``
* ``structural_temperature``

Displacement
------------
Displacement is the DOF solution for structural analyses.  Its
:class:`Displacement <ansys.dpf.post.displacement.Displacement>` object can be
with the following:

.. code:: python

    # Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the displacement result object

    >>> displacement = solution.displacement()

Being a DOF solution, its location argument **must be nodal** (inferred above).

The displacement ``Result`` corresponds to a vector field. To obtain the scalar
components of this field, e.g., y-components, access the following subresult:

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

For more details, see :ref:`ref_api_result_data`.


Stress
------
The :class:`Stress <ansys.dpf.post.stress.Stress>` result object can be
with the following:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    # Instantiate the stress result object

    >>> stress = solution.stress()

The stress result corresponds to a tensor field. To obtain the scalar
components of this field, such as the normal y-stresses, access the following
subresult:

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

Other components, as well as whole tensor data, can be queried accordingly.
For more details, see :ref:`ref_api_result_data`.


Strain (elastic, plastic)
-------------------------
The :class:`ElasticStrain <ansys.dpf.post.strain.ElasticStrain>` and
:class:`PlasticStrain <ansys.dpf.post.strain.PlasticStrain>` result objects can
be obtained with:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the strain result objects

    >>> elastic_strain = solution.elastic_strain()
    >>> plastic_strain = solution.plastic_strain()

A strain result corresponds to a tensor field. To obtain the scalar components
of this field, such as the shear xy-strains, access the following subresult:

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

Other components, as well as whole tensor data, can be queried accordingly.
For more details, see :ref:`ref_api_result_data`.


Structural temperature
----------------------
The :class:`StructuralTemperature
<ansys.dpf.post.temperature.StructuralTemperature>` result object can be
obtained with the following:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate the structural temperature result object

    >>> structural_temperature = solution.structural_temperature()

To access the temperature scalar field use the following:

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
Other miscellaneous :class:`ansys.dpf.post.misc_results.MecanicMisc` result
objects may be available in the ``Solution``. For example:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Get the result data

    >>> acceleration = solution.misc.nodal_acceleration()

All **keyword arguments** are available for miscellaneous results (see
:ref:`ref_result_keywords`), except location. Notably, some subresults
may be available as keyword arguments, for example, the scalar
components of acceleration:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Get the result data

    >>> acceleration = solution.misc.nodal_acceleration(subresult="Y")

Verify the result file contents to determine available queries (see
:ref:`user_guide_accessing_file_metadata`).


Thermal/Electric (``*.rth``) result files
=========================================
Upon loading a ``Solution`` object from a thermal/electric analysis
results file (\*.rth), the following ``Result`` objects may be
queried:

* ``temperature``
* ``heat_flux``
* ``electric_field``
* ``electric_potential``

Temperature
-----------
Temperature is the DOF solution for thermal analyses. The :class:`Temperature
<ansys.dpf.post.temperature.Temperature>` result object can be obtained as
follows:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.steady_therm)

    Instantiate the temperature result object

    >>> temperature = solution.temperature()

Being a DOF solution, its location argument **must be nodal**
(inferred above).  The scalar field can be obtained directly:

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
The :class:`HeatFlux <ansys.dpf.post.temperature.HeatFlux>` result object can
be obtained with the following:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.steady_therm)

    Instantiate the heat_flux result object

    >>> heat_flux = solution.heat_flux()

The heat flux result corresponds to a vector field. To obtain the scalar
components of this field, e.g., x-components, access the following subresult:

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

Other components can be queried accordingly.
For more details, see :ref:`ref_api_result_data`.


Electric field
--------------
The :class:`ElectricField <ansys.dpf.post.electric_results.ElectricField>`
result object can be obtained with the following:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.electric_therm)

    Instantiate the electric field result object

    >>> electric_field = solution.electric_field()

The electric field result corresponds to a vector field. To
obtain the scalar components of this field, such as the x-components, access
the following subresult:

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

For more details, see :ref:`ref_api_result_data`.


Electric potential
------------------
The :class:`ElectricPotential
<ansys.dpf.post.electric_results.ElectricPotential>` result object can be
with:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.steady_therm)

    Instantiate the electric potential result object

    >>> electric_potential = solution.electric_potential()

It represents a scalar field whose values can be obtained with:

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
