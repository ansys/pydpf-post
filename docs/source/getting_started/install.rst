.. _installation:

************
Installation
************

Install using ``pip``
---------------------

The standard package installer for Python is `pip <https://pypi.org/project/pip/>`_.

To use PyDPF-Post with Ansys 2021 R1 or later, install the latest version
with this command:

.. code::

   pip install ansys-dpf-post

To install PyDPF-Post with its optional plotting functionalities, use:

.. code::
   pip install ansys-dpf-post[plotting]


Install without internet
------------------------

If you are unable to install PyDPF-Post on the host machine using ``pip`` due to
network isolation, download the wheelhouse corresponding to your platform and Python interpreter version
for the latest release of PyDPF-Post from the assets section of the `latest PyDPF-Post release on GitHub <https://github.com/pyansys/pydpf-post/releases/latest>`_. 

The wheelhouse is a ZIP file containing Python wheels for all the packages PyDPF-Post requires to run.
To install PyDPF-Post using the downloaded wheelhouse, unzip the wheelhouse to a local directory,
then use the following command from within this local directory:
.. code::

   pip install --no-index --find-links=. ansys-dpf-post

Beware that PyDPF-Post wheelhouses do not include the optional plotting dependencies.
To allow for plotting capabilities, also download the wheel corresponding to your platform and Python interpreter version
for `PyVista <https://pypi.org/project/pyvista/#files>`_, then place it in the same previous local directory and run the command above.

Install in development mode
---------------------------

If you want to edit and potentially contribute to DPF-Post,
clone the repository and install it using ``pip`` with the ``-e``
development flag:

.. include:: ../pydpf-post_clone_install.rst
