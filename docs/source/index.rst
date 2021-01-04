================
PyANSYS DPF-Post
================

The Data Processing Framework (DPF) is designed to provide numerical
simulation users/engineers with a toolbox for accessing and
transforming simulation data. DPF can access data from solver result
files as well as several neutral formats (csv, hdf5, vtk,
etc.). Various operators are available allowing the manipulation and
the transformation of this data.

The Python ``ansys.dpf.post`` module provides an easy Python interface 
enabling rapid post-processing work in a simple way, without ever leaving a
Python environment. It is based on a solution instance and result objects. 

To learn more about PDF mechanisms, please try to use the 
``ansys.dpf.core`` module. 


Brief Demo
~~~~~~~~~~
Opening a result file generated from Ansys workbench or MAPDL is as easy as:

.. code:: python

    >>> from ansys.dpf import post
    >>> solution = post.load_solution('../../tests/testfiles/model_with_ns.rst')
	>>> displacement = solution.displacement(node_scoping = 3)
	>>> result = displacement.vector

See the :ref:`gallery` for detailed examples.


Key Features
~~~~~~~~~~~~

**Computation Efficiency**

As the DPF-Post module is based on DPF Framework (that been modernly developed by taking advantages of new hardware architectures), 
it relies on its computation efficiency.

**Easy to use**

The API of the Post module has been developped in order to make easy post-processing steps easier. It is also an introduction of 
the DPF's concepts. 


.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   getting_started/index
   user_guide
   api/index
   examples/index
   contributing
