.. _user_guide_extending_to_core:

**********************************************
Extending DPF-Post capabilities using DPF-Core
**********************************************

The DPF-Post module provides an simplified Python interface to DPF, 
thus enabling rapid post-processing. Another DPF module can be used 
with a Python interface: DPF-Core. The Post module is based on the 
concepts of the Core one. 

The Data Processing Framework (DPF-Core) can access data from solver result 
files as well as several neutral formats. Various **operators** are available allowing 
the manipulation and the transformation of this data. It allows simple and/or 
complex evaluations by chaining operators. The data in DPF is defined based 
on physics agnostic mathematical quantities described in a self-sufficient 
entity called **field**. The **fields container** contain field(s).

This module is based on the DPF-Core module, that allows more 
capabilities. A few will be introduced below.


Export VTK
----------

The following code shows how to proceed to **export a fields container in VTK format**:

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
	>>> # the result fields_container is extracted from the result data
	>>> fields_container = norm.result_fields_container
	
	>>> # Now the Core API need to be imported
	>>> from ansys.dpf import core
	>>> # a dedicated operator need to be instantiated
	>>> vtk_operator = core.Operator("vtk_export")
	>>> # connections must be set
	>>> vtk_operator.inputs.mesh.connect(solution.mesh)
	>>> vtk_operator.inputs.file_path.connect("vtk_example.vtk")
	>>> vtk_operator.inputs.fields1.connect(fields_container)
	>>> # the operator need to be run
	>>> vtk_operator.run()


Export HDF5
-----------

The following code shows how to proceed to **export a fields container in HDF5 format**:

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
	>>> # the result fields_container is extracted from the result data
	>>> fields_container = norm.result_fields_container
	
	>>> # Now the Core API need to be imported
	>>> from ansys.dpf import core
	>>> # a dedicated operator need to be instantiated
	>>> h5_operator = core.Operator("serialize_to_hdf5")
	>>> # connections must be set
	>>> h5_operator.inputs.mesh.connect(solution.mesh)
	>>> h5_operator.inputs.file_path.connect("hdf5_example.h5")
	>>> h5_operator.inputs.data.connect(fields_container)
	>>> # the operator need to be evaluated
	>>> h5_operator.eval()

