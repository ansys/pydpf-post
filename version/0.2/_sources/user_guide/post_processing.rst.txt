.. _user_guide_post_processing:

********************
Load the result file
********************

The :class:`DpfSolution <ansys.dpf.post.dpf_solution.DpfSolution>` object is
a central element of DPF-Post. This object is the entry point for browsing
the contents of a result file.

**On Windows**

You can load the result file with:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution('C:/Users/user/file.rst')


**On Linux**
You can load the result file with:
    
.. code:: python

    >>> solution = post.load_solution('/home/user/file.rst')


For a more detailed example on interacting with the
:class:`DpfSolution <ansys.dpf.post.dpf_solution.DpfSolution>` object,
see :ref:`ref_basics`.
