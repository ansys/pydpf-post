.. _user_guide_post_processing:

************************
Loading the Results File
************************

The DPF-Post package provides a Python interface to DPF (Data
Processing Framework) that is streamlined for post-processing.

A central concept in DPF-Post is the :class:`DpfSolution
<ansys.dpf.post.dpf_solution.DpfSolution>` class. It is the entry point for
browsing the contents of a result file.  It can be obtained as follows:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution('C:/Users/user/file.rst')

    # Or on linux

    >>> solution = post.load_solution('/home/user/file.rst')

For a a more detailed example on interacting with the :class:`DpfSolution
<ansys.dpf.post.dpf_solution.DpfSolution>` class, see :ref:`ref_basics`.
