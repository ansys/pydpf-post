.. _user_guide_troubleshooting:

===============
Troubleshooting
===============
This section explains how to resolve the most common issues encountered with ``pydpf-post``.
It also includes suggestions for improving scripts.

Using the Solution
------------------

Invalid UTF-8 warning/issue 
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Trying to load a solution with: 

.. code-block:: default

    from ansys.dpf import post
    solution = post.load_solution('file.rst')

The following error can be raised: 

.. code-block:: default

    [libprotobuf ERROR C:\.conan\897de8\1\protobuf\src\google\protobuf\wire_format_lite.cc:578] 
    String field 'ansys.api.dpf.result_info.v0.ResultInfoResponse.user_name' contains invalid UTF-8 
    data when serializing a protocol buffer. Use the 'bytes' type if you intend to send raw bytes.

This will prevent the solution to be accessed. To avoid a such inconvenience, please ensure to work with 
a version higher than 0.2.1 for the ansys-dpf-post module, and higher than 0.3.2 of the ansys-dpf-core module.
In this case, a warning will still be raised, but it should not prevent to use the solution anymore. 

The solution needs to set a physics_type and an analysis type to allow the result file to be read. This
can be specified using the "physics_type" and "analysis_type" arguments of the load_solution method: 

.. code-block:: default

    from ansys.dpf import post
    solution = post.load_solution('file.rst', physics_type='mecanic', analysis_type='transient')