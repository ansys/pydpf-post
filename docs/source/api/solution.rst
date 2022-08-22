.. _ref_api_solution:

***************
``DpfSolution``
***************

The :class:`DpfSolution <ansys.dpf.post.dpf_solution.DpfSolution>` class is
the entry point for browsing the contents of the results file. Use the
following code to load the file and obtain the instance.

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

.. autosummary::
   :toctree: _autosummary

   ansys.dpf.post.load_solution
   ansys.dpf.post.dpf_solution.DpfSolution
