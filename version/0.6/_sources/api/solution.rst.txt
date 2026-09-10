.. _ref_api_solution:

*********************
``DpfSolution`` class
*********************

The :class:`DpfSolution <ansys.dpf.post.dpf_solution.DpfSolution>` class is
the entry point for browsing the contents of the result file.

This code shows how to load a result file:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

.. autosummary::
   :toctree: _autosummary

   ansys.dpf.post.load_solution
   ansys.dpf.post.dpf_solution.DpfSolution
