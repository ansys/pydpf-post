.. _installation:

************
Installation
************

Install using ``pip``
---------------------

The standard package installer for Python is `pip <https://pypi.org/project/pip/>`_.

To use PyDPF-Post with Ansys 2021 R1 or later, install the latest version
with this command:

.. code:: bash

   pip install ansys-dpf-post

PyDPF-Post plotting capabilities require an installation of `PyVista <https://pyvista.org/>`_.
To install PyDPF-Post with its optional plotting functionalities, use this command:

.. code:: bash

   pip install ansys-dpf-post[plotting]

For more information about PyDPF-Post plotting capabilities, see :ref:`user_guide_plotting`.


Install without internet
------------------------

If you are unable to install PyDPF-Post on the host machine using ``pip`` due to
network isolation, you can download the wheelhouse corresponding to your platform
and Python interpreter version. To obtain the latest release, go to the **Assets** section
for the `latest PyDPF-Post release <https://github.com/ansys/pydpf-post/releases/latest>`_ on GitHub.

The wheelhouse is a ZIP file containing Python wheels for all the packages PyDPF-Post requires to run.
To install PyDPF-Post using the downloaded wheelhouse, unzip the wheelhouse to a local directory and
then run the following command from within this local directory:

.. code:: bash

   pip install --no-index --find-links=. ansys-dpf-post

Beware that PyDPF-Post wheelhouses do not include the optional plotting dependencies.
To allow for plotting capabilities, also download the `PyVista wheel <https://pypi.org/project/pyvista/#files>`_
and unzip it to the same local directory before running the preceding command again.


Install in development mode
---------------------------

If you want to edit and potentially contribute to PyDPF-Post,
clone the repository and install it using ``pip`` with the ``-e``
development flag:

.. include:: ../pydpf-post_clone_install.rst


Check the installation
----------------------

Run the following Python code to verify your PyDPF-Post installation:

.. code:: bash

   from ansys.dpf import post
   from ansys.dpf.post import examples
   simulation = post.load_simulation(examples.simple_bar)
   print(simulation)
