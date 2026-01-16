.. _user_guide_troubleshooting:

===============
Troubleshooting
===============
This section explains how to resolve the most common issues encountered with PyDPF-Post.
It also includes suggestions for improving scripts.

Invalid UTF-8 error
~~~~~~~~~~~~~~~~~~~
Assume that you are using this code to load a solution: 

.. code-block:: default

    from ansys.dpf import post
    solution = post.load_solution('file.rst')

This error might be raised: 

.. code-block:: default

    [libprotobuf ERROR C:\.conan\897de8\1\protobuf\src\google\protobuf\wire_format_lite.cc:578] 
    String field 'ansys.api.dpf.result_info.v0.ResultInfoResponse.user_name' contains invalid UTF-8 
    data when serializing a protocol buffer. Use the 'bytes' type if you intend to send raw bytes.

This error prevents the solution from being accessed. To avoid this error, ensure that you are using
a PyDPF-Post version later than 0.2.1 and a PyDPF-Core version later than 0.3.2.
In this case, a warning might still be raised, but it should not prevent the solution from being accessed.

For the result file to be read, you must set the ``physics_type`` and ``analysis type`` arguments for the
``load_solution ()`` method:

.. code-block:: default

    from ansys.dpf import post
    solution = post.load_solution('file.rst', physics_type='mecanic', analysis_type='transient')
