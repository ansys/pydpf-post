.. _ref_api_solution:

***************
Solution object
***************

The ``Solution`` object is the entry point for browsing the contents
of the results file.  Use the following code to load the file and
obtain the instance.

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)
	

.. autoclass:: ansys.dpf.post.dpf_solution.DpfSolution
    :members:

