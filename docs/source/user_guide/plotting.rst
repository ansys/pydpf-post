.. _user_guide_plotting:

********
Plotting
********
The DPF-Post package also provides functionality to conveniently plot your result.
It suffices to load the ``Solution`` object with the results file, request a ``Result`` object, 
and obtain the scalar field of interest.  The ``plot_contour`` method will be ready to render it.


Total deformation 
-----------------

Use the following code to **plot the total deformation** 
(norm of the displacement vector field):

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

Use the following code to **plot the normal x-component of stress**:

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

