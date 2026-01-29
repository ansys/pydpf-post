.. _ref_api_result_object:

****************
``Result`` class
****************

The :class:`Result <ansys.dpf.post.result_object.Result>` class provides access
to different types of result data.

This code shows how you can access result data for displacements, stresses, and
elastic strains:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Displacement result object

    >>> displacement = solution.displacement()

    Stress result object

    >>> stress = solution.stress()

    Elastic strain result object

    >>> elastic_strain = solution.elastic_strain()

DPF-Post features separate ``Result`` classes for each supported result type. For information on
the available types and their interfaces, see :ref:`user_guide_accessing_results`.

.. currentmodule:: ansys.dpf.post

.. autosummary::
   :toctree: _autosummary

   result_object.Result
   displacement.Displacement
   electric_results.ElectricField
   electric_results.ElectricPotential
   misc_results.MecanicMisc
   strain.ElasticStrain
   strain.PlasticStrain
   stress.Stress
   temperature.HeatFlux
   temperature.StructuralTemperature
   temperature.Temperature
