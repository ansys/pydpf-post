.. _user_guide_accessing_file_metadata:

************************
Browsing result metadata
************************

The ``Solution`` object is the entry point for browsing the contents of a 
result file (see :ref:`user_guide_post_processing`). It provides access 
to important metadata such as the analysis type and listings of available results. 

This information is accessible as follows.

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # access the result information 
	>>> result_info = solution.get_result_info()
	>>> print(result_info)
	
Four different analysis types are presently supported in DPF-Post:

* static analysis
* modal analysis 
* harmonic analysis
* transient analysis 
	
Mesh
----
The mesh can also be accessed from the ``Solution`` object. 

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # get the mesh
	>>> mesh = solution.mesh
	>>> print(mesh)

