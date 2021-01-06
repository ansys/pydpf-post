.. _ref_api_index:

API Reference
=============
Details of the DPF API.


Solution Object
---------------
The solution object instantiates an object that is built on the result
file.  Use the following code to instantiate a solution object.

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)


Result Object
-------------
The result object can be manipulated to get different result data.

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


ResultData Class
----------------
The ``ResultData`` class is instantiated from a result object.  It
enables easy access to the data contained within the result
object. The following example shows how to get an instance of
``ResultData``.
    
.. code:: python

    >>> from ansys.dpf import post
    >>> solution = post.load_solution(r'../../../testfiles/tests/model_with_ns.rst')
    >>> displacement = solution.displacement()
    >>> result = displacement.vector
    >>> # access the data
    >>> vector.get_data_at_field(0)
    
.. autoclass:: ansys.dpf.post.result_data.ResultData
    :members:
