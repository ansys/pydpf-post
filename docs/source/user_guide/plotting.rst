.. _user_guide_plotting:

********
Plotting
********
The DPF-Post package also provides functionality to conveniently plot
your result.  It suffices to load the ``Solution`` object with the
results file, request a ``Result`` object, and obtain the scalar field
of interest.  The ``plot_contour`` method will be ready to render it.


Total deformation 
-----------------

Use the following code to **plot the total deformation** 
(norm of the displacement vector field):

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate a displacement result object 

    >>> displacement = solution.displacement()
    >>> norm = displacement.norm  # this is the result data (data container)

    Plot the result data

    >>> norm.plot_contour()


Normal Stresses
---------------

Use the following code to **plot the normal x-component of stress**:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate a stress result object

    >>> stress = solution.stress()
    >>> s_xx = stress.xx  # this is the result data (data container)

    Plot the result data

    >>> s_xx.plot_contour()

