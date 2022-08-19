.. _user_guide_extending_to_core:

*****************************************
Using DPF-Core for More General Operators
*****************************************

The DPF-Post package provides a Python interface to DPF (Data
Processing Framework) that is streamlined for post-processing. It also
works hand-in-hand with the wider DPF-Core package that is especially
suited for representing more general data transformations with the use
of powerful, scalable operators.

DPF-Core is equipped to access data from Ansys solver result files as well as
third-party formats. Various `operators
<https://dpfdocs.pyansys.com/operator_reference.html>`_ are available to
facilitate data manipulation and transformation. DPF-Core allows arbitrarily
complex workflows to be assembled from simpler building blocks, by chaining
together operators. The data in DPF is represented by physics-agnostic
mathematical quantities described in a self-sufficient entity called a
:class:`Field <ansys.dpf.core.field.Field>`.

Because DPF-Post is based on DPF-Core, it is helpful to see their
interplay to appreciate their range of capabilities, as done in the
following examples.


Export VTK
----------

The following code shows how to **export a fields container in VTK format**:

.. code:: python

    Instantiate the solution object.

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate a result object. This is a displacement result.

    >>> displacement = solution.displacement()

    This is the result data (data container).

    >>> norm = displacement.norm

    The result fields_container is extracted from the result data.

    >>> fields_container = norm.result_fields_container
    
    Now the Core API need to be imported.

    >>> from ansys.dpf import core

    A dedicated operator needs to be instantiated.

    >>> vtk_operator = core.Operator("vtk_export")

    Connections must be set.

    >>> vtk_operator.inputs.mesh.connect(solution.mesh)
    >>> vtk_operator.inputs.file_path.connect("vtk_example.vtk")
    >>> vtk_operator.inputs.fields1.connect(fields_container)

    Run the operator.

    >>> vtk_operator.run()


Export HDF5
-----------

The following code shows how to proceed to **export a fields container
in HDF5 format**:

.. code:: python

    Instantiate the solution object 

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Instantiate a result object. This is a displacement result.

    >>> displacement = solution.displacement()

    This is the result data (data container)

    >>> norm = displacement.norm

    Extract the result fields_container from the result data.

    >>> fields_container = norm.result_fields_container
    
    Now the Core API needs to be imported.

    >>> from ansys.dpf import core

    Initialize a dedicated operator.

    >>> h5_operator = core.Operator("serialize_to_hdf5")

    Set the connection.

    >>> h5_operator.inputs.mesh.connect(solution.mesh)
    >>> h5_operator.inputs.file_path.connect("hdf5_example.h5")
    >>> h5_operator.inputs.data.connect(fields_container)

    Evaluate the the operator.

    >>> h5_operator.eval()

