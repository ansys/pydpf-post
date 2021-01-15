.. _ref_api_result_data:

****************
ResultData Class
****************

The ``ResultData`` class is instantiated from a result object.  It
enables easy access to the data contained within the result
object. The following example shows how to get an instance of
``ResultData``.
    
.. code:: python

    >>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
    >>> displacement = solution.displacement()
    >>> result = displacement.vector
    >>> # access the data
    >>> vector.get_data_at_field(0)
    
.. autoclass:: ansys.dpf.post.result_data.ResultData
    :members:
