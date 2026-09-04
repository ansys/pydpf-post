.. _user_guide_post_processing:

********************
Load the result file
********************

The :class:`Simulation <ansys.dpf.post.simulation.Simulation>` object is
a central element of PyDPF-Post. This object is the entry point for browsing
the contents of a result file.

**On Windows**

You can load the result file with:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> simulation = post.load_simulation('C:/Users/user/file.rst')


**On Linux**
You can load the result file with:
    
.. code:: python

    >>> simulation = post.load_simulation('/home/user/file.rst')


For a more detailed example on interacting with the
:class:`Simulation <ansys.dpf.post.simulation.Simulation>` object,
see :ref:`ref_basics`.
