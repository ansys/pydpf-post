.. _ref_api_result_object:

*************
Result object
*************

The ``result object`` can be manipulated to get different result data.

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)
    >>> # Displacement result object
    >>> displacement = solution.displacement()
    >>> # Stress result object
    >>> stress = solution.stress()
    >>> # Elastic strain result object
    >>> elastic_strain = solution.elastic_strain()

A result object is available for each main DPF-Post result type. 
See the following list to know which result can be accessed this way. 

List of available result objects using DPF-Post API:
- displacement
- stress
- elastic_strain
- plastic_strain
- structural_temperature

See the :ref:`_user_guide_accessing_results` section to get more information 
about the API of each specific result object. 