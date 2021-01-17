.. _user_guide_post_processing:

************************
Loading the Results File
************************

The DPF-Post package provides a Python interface to DPF (Data
Processing Framework) that is streamlined for post-processing.

A central concept in DPF-Post is the ``Solution`` object. It is the
entry point for browsing the contents of a result file.  It can be
obtained as follows:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution('C:/Users/user/file.rst')

    or on linux

    >>> model = dpf.Model('/home/user/file.rst')

For a a more detailed example on interacting with the ``Solution``
object, see :ref:`ref_basics`.
