.. _ref_api_solution:

**************
Solution Class
**************

The ``solution`` object instantiates an object that is built on the result
file.  Use the following code to instantiate a solution object.

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)
	

.. autoclass:: ansys.dpf.post.dpf_solution.DpfSolution
    :members:

