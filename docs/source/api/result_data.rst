.. _ref_api_result_data:

*****************
ResultData object
*****************

The ``ResultData`` object is obtained from a ``Result``.  It provides
access to the result values contained within, as shown in the
following example:
    
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
