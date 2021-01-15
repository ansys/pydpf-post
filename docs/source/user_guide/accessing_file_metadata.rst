.. _user_guide_accessing_file_metadata:

*******************************
Accessing results file metadata
*******************************

The DPF-Post module provides an simplified Python interface to DPF, 
thus enabling rapid post-processing. To proceed, a solution object 
must be instantiated first.
Refer to the previous part to proceed (:ref:`user_guide_post_processing`). 

Once the solution is instantiated, it is possible to read the file metadata. 

Analysis type
-------------
The result information will provide the result file's analysis type.

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # access the result information 
	>>> result_info = solution.get_result_info()
	>>> print(result_info)
	
There are four different analysis types covered using the DPF-Post module:
- static analysis
- modal analysis 
- harmonic analysis
- transient analysis 

Available results
-----------------
The result information will also provide the available results. 

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # access the result information 
	>>> result_info = solution.get_result_info()
	>>> print(result_info)
	
Mesh
----
The mesh can be accessed from the solution. 

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # get the mesh
	>>> mesh = solution.mesh
	>>> print(mesh)

