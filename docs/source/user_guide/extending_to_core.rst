.. _user_guide_extending_to_core:

**********************
Use DPF-Core operators
**********************

DPF-Post is based on DPF-Core, but it is streamlined for postprocessing.
Because DPF-Post works hand-in-hand with DPF-Core, you can use the powerful,
scalable  `operators <https://dpfdocs.pyansys.com/operator_reference.html>`_
in DPF-Core to facilitate data manipulation and more general data transformations.

DPF-Core can access data from Ansys solver result files as well as from
third-party file formats. With DPF-Core, you can assemble complex workflows
from simpler building blocks by chaining them together with operators.
The data in DPF is represented by physics-agnostic mathematical quantities
described in a self-sufficient entity called a :class:`Field <ansys.dpf.core.field.Field>`.

To show the range of DPF-Core and DPF-Post capabilities, the following
examples show they work together.


Export to VTK file
------------------

This code shows how to export a fields container to a VTK file:

.. code:: python

    # Instantiate the solution object.

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate a result object. This is a displacement result.

    >>> displacement = solution.displacement()

    This is the result data (data container).

    >>> norm = displacement.norm

    Extract the fields_container result from the result data.

    >>> fields_container = norm.result_fields_container
    
    Import the DFP-Core API.

    >>> from ansys.dpf import core

    # Instantiate a dedicated operator

    >>> vtk_operator = core.Operator("vtk_export")

    Set connections.

    >>> vtk_operator.inputs.mesh.connect(solution.mesh)
    >>> vtk_operator.inputs.file_path.connect("vtk_example.vtk")
    >>> vtk_operator.inputs.fields1.connect(fields_container)

    Run the operator.

    >>> vtk_operator.run()


Export to HDF5 file
-------------------

This code shows how to export the same fields container
to an HDF5 file:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate a result object. This is a displacement result.

    >>> displacement = solution.displacement()

    This is the result data (data container).

    >>> norm = displacement.norm

    Extract the fields_container result from the result data.

    >>> fields_container = norm.result_fields_container
    
    Import the DPF-Core API.

    >>> from ansys.dpf import core

    Initialize a dedicated operator.

    >>> h5_operator = core.Operator("serialize_to_hdf5")

    Set the connection.

    >>> h5_operator.inputs.mesh.connect(solution.mesh)
    >>> h5_operator.inputs.file_path.connect("hdf5_example.h5")
    >>> h5_operator.inputs.data.connect(fields_container)

    Evaluate the operator.

    >>> h5_operator.eval()

