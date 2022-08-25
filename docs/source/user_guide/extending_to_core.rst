.. _user_guide_extending_to_core:

***************************************
Use DPF-Core for more general operators
***************************************

DPF-Post, streamlined for postprocessing, works hand-in-hand
with DPF-Core. This broader package provides powerful, scalable
operators for more general data transformations.

DPF-Core can access data from Ansys solver result files as well as
from third-party file formats. It provides various `operators
<https://dpfdocs.pyansys.com/operator_reference.html>`_ that
facilitate data manipulation and transformation. With DPF-Core, you
can assemble complex workflows from simpler building blocks by chaining
them together with operators. The data in DPF is represented by physics-agnostic
mathematical quantities described in a self-sufficient entity called a
:class:`Field <ansys.dpf.core.field.Field>`.

DPF-Post is based on DPF-Core. Thus, to appreciate the range of their
capabilities, the following examples show how the effect they have on each
other.


Export VTK
----------

This code shows how to export a fields container in VTK format:

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


Export HDF5
-----------

This code shows how to subsequently export the fields container
in HDF5 format:

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

