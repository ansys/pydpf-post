.. _user_guide_troubleshooting:

===============
Troubleshooting
===============
This section explains how to resolve the most common issues encountered with PyDPF-Post.
It also includes suggestions for improving scripts.


Known Issues List
~~~~~~~~~~~~~~~~~
Please refer to :ref:`Known Issues List <ref_kil>` for a list of known issues for each version of DPF.


Installation
~~~~~~~~~~~~
When pip installing older versions of the PyDPF libraries, an error might occur stating:

.. code-block:: default

    ``"'python_requires' must be a string containing valid version specifiers;
    Invalid specifier: '>=3.7.*'"``.

In this case, modify your Python environment to use a version of the ``setuptools`` library strictly
older than 67.0.0 using the command below:

.. code-block:: default

    pip uninstall -y setuptools; pip install setuptools<67.0.0


Auto-completion
~~~~~~~~~~~~~~~
Depending on your scripting environment, auto-completion might not work correctly when using the
``load_simulation()`` method. This method is intended as a helper which can detect automatically
the physics type and analysis type. To get auto-completion to work in all circumstances, instantiate
 the right :ref:`Simulation <ansys.dpf.post.simulation.Simulation` sub-class directly using its
constructor:

.. code-block:: default

    from ansys.dpf import post

    static_mechanical_simulation = post.StaticMechanicalSimulation('file.rst')
    # or
    transient_mechanical_simulation = post.TransientMechanicalSimulation('file.rst')
    # or
    modal_mechanical_simulation = post.ModalMechanicalSimulation('file.rst')
    # or
    harmonic_mechanical_simulation = post.HarmonicMechanicalSimulation('file.rst')


Invalid UTF-8 error
~~~~~~~~~~~~~~~~~~~
Assume that you are using this code to load a simulation result:

.. code-block:: default

    from ansys.dpf import post
    simulation = post.load_simulation('file.rst')

This error might be raised: 

.. code-block:: default

    [libprotobuf ERROR C:\.conan\897de8\1\protobuf\src\google\protobuf\wire_format_lite.cc:578] 
    String field 'ansys.api.dpf.result_info.v0.ResultInfoResponse.user_name' contains invalid UTF-8 
    data when serializing a protocol buffer. Use the 'bytes' type if you intend to send raw bytes.

This error prevents the simulation from being accessed. To avoid this error, ensure that you are using
a PyDPF-Post version later than 0.2.1 and a PyDPF-Core version later than 0.3.2.
In this case, a warning might still be raised, but it should not prevent the simulation from being accessed.

For the result file to be read, you must set the ``physics_type`` and ``analysis type`` arguments for the
``load_solution()`` method:

.. code-block:: default

    from ansys.dpf import post
    solution = post.load_solution('file.rst', physics_type='mechanical', analysis_type='transient')
