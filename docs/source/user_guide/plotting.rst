.. _user_guide_plotting:

********
Plotting
********

The DPF-Post module provides an simplified Python interface to DPF, 
thus enabling rapid post-processing. To proceed, a solution object 
must be instantiated first.
Refer to the part to proceed (:ref:`user_guide_post_processing'). 

Once the solution is instantiated, it is possible to access to the 
contained results.
Refer to the previous part to proceed (:ref:`user_guide_accessing_results'). 

Results can be manipulated as result data objects. The API allows to 
plot the wanted result on its mesh. 


Total deformation 
-----------------

The following code shows how to proceed to **plot the total deformation** 
(norm of the displacement vector):

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate a result object 
	>>> # here is called a displacement result
	>>> displacement = solution.displacement()
	>>> # here is called the result data (data container)
	>>> norm = displacement.norm
	>>> # here the result data is plotted
	>>> norm.plot_contour()


Stress XX component
-------------------

The following code shows how to proceed to **plot the stress XX component**:

.. code:: python

	>>> # instantiate the solution object 
	>>> from ansys.dpf import post
	>>> from ansys.dpf.post import examples
	>>> solution = post.load_solution(examples.multishells_rst)
	>>> # instantiate a result object 
	>>> # here is called a stress result
	>>> stress = solution.stress()
	>>> # here is called the result data (data container)
	>>> s_xx = stress.xx
	>>> # here the result data is plotted
	>>> s_xx.plot_contour()

