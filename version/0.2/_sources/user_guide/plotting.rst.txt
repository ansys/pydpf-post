.. _user_guide_plotting:

************
Plot results
************

DPF-Post provides functionality for plotting results. Here is a summary of
the steps:

#. Load the :class:`DpfSolution <ansys.dpf.post.dpf_solution.DpfSolution>` object
   with the result file.
#. Request a :class:`Result <ansys.dpf.post.result_object.Result>` object and
   obtain the scalar field of interest.
#. Use the :func:`plot_contour <ansys.dpf.post.result_data.ResultData.plot_contour>`
   method to render it.

Subsequent sections provide some plotting examples.

Total deformation 
-----------------

You can plot the total deformation (norm of the displacement vector field) with:

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


Normal stresses
---------------

You can plot the normal x-component of stress with:

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

