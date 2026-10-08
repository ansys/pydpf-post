.. _installation:

************
Installation
************
After all :ref:`software requirements <software_requirements>` are installed, you can install PyDPF-Post.

Install using ``pip``
---------------------

The standard package installer for Python is `pip <https://pypi.org/project/pip/>`_.

To use PyDPF-Post with Ansys 2021 R1 or later, install the latest version
with this command:

.. code:: bash

   pip install ansys-dpf-post

Alternatively, you can also clone and install the latest version of PyDPF-Post from its GitHub
repository with these commands:

.. code:: bash

   git clone https://github.com/ansys/pydpf-post
   cd pydpf-post
   pip install . --userpip install ansys-dpf-post

PyDPF-Post plotting capabilities require an installation of `PyVista <https://pyvista.org/>`_.
To install PyDPF-Post with its optional plotting functionalities, use this command:

.. code:: bash

   pip install ansys-dpf-post[graphics]

.. warning::

   ``pip install ansys-dpf-post[plotting]`` is equivalent to the previous command, however, the "plotting" target
   only remains valid for legacy reasons and will soon be deprecated. Users are encouraged to use the "graphics"
   target instead.

For more information about PyDPF-Post plotting capabilities, see :ref:`user_guide_plotting`.


Install without internet
------------------------

If you are unable to install PyDPF-Post on the host machine using ``pip`` due to
network isolation, you can download the wheel file corresponding to your platform
and Python interpreter version. To obtain the latest release, go to the **Assets** section
for the `latest PyDPF-Post release <https://github.com/ansys/pydpf-post/releases/latest>`_ on GitHub.

A Python wheel file is essentially a ZIP archive containing Python wheel files for all the packages
that PyDPF-Post requires to run. To install PyDPF-Post using the downloaded wheel file, unzip this file
to a local directory and then run this command from within this local directory:

.. code:: bash

   pip install --no-index --find-links=. ansys-dpf-post

Beware that PyDPF-Post wheel files do not include the optional plotting dependencies.
To allow for plotting capabilities, also download the `PyVista wheel file <https://pypi.org/project/pyvista/#files>`_,
unzip it to the same local directory, and run the preceding command again.

Check the installation
----------------------

Run this Python code to verify your PyDPF-Post installation:

.. code:: bash

   from ansys.dpf import post
   from ansys.dpf.post import examples
   simulation = post.load_simulation(examples.simple_bar)
   print(simulation)

