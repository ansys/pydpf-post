.. _ref_api_result_object:

*******
Results
*******

The :class:`Result <ansys.dpf.post.result_object.Result>` class provides access
to different types of result data.  See the example below for access to
displacements, stresses and elastic strains.

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

DPF-Post features separate result classes for each supported result type. See
the :ref:`user_guide_accessing_results` section for information on the
available types and their interfaces.

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
