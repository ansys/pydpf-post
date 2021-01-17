.. _ref_api_result_object:

*************
Result Object
*************

The ``Result`` object provides access to different types of result
data.  See the example below for access to displacements, stresses and
elastic strains.

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

DPF-Post features separate ``Result`` classes for each supported
result type. See the :ref:`user_guide_accessing_results` section for
information on the available types and their interfaces.
