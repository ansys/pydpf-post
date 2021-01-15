.. _user_guide_post_processing:

****************************
Post-processing results file
****************************

The DPF-Post module provides an simplified Python interface to DPF, 
thus enabling rapid post-processing.

To proceed, a solution object must be instantiated first.
It is an object that is built on the result file. 
Use the following code to instantiate a solution object.

To create a ``Solution`` instance, import ``dpf-post`` and load a file.  The
path provided must be an absolute path or a path relative to the file you are 
writting.

.. code:: python

	from ansys.dpf import post
	from ansys.dpf.post import examples
	solution = post.load_solution('C:/Users/user/file.rst')
	# or on linux
	model = dpf.Model('/home/user/file.rst')


For a full and easy example using the solution, see :ref:`_ref_basics`.