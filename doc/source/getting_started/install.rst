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

Or, you can clone and install the latest version of PyDPF-Post from its GitHub
repository with these commands:

.. code:: bash

   git clone https://github.com/ansys/pydpf-post
   cd pydpf-post
   pip install . --userpip install ansys-dpf-post

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

The wheelhouse is a ZIP file containing Python wheels for all the packages that PyDPF-Post requires to run.
To install PyDPF-Post using the downloaded wheelhouse, unzip the wheelhouse to a local directory and
then run this command from within this local directory:

.. code:: bash

   pip install --no-index --find-links=. ansys-dpf-post

Beware that PyDPF-Post wheelhouses do not include the optional plotting dependencies.
To allow for plotting capabilities, also download the `PyVista wheel <https://pypi.org/project/pyvista/#files>`_,
unzip it to the same local directory, and run the preceding command again.

Check the installation
----------------------

Run this Python code to verify your PyDPF-Post installation:

.. code:: bash

   from ansys.dpf import post
   from ansys.dpf.post import examples
   simulation = post.load_simulation(examples.simple_bar)
   print(simulation)

Brief demo
----------

Provided you have Ansys 2023 R1 or later installed, a DPF server automatically starts
once you start using PyDPF-Post.

To load a simulation for a MAPDL result file to extract and
postprocess results, use this code:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> simulation = post.load_simulation(examples.download_crankshaft())
    >>> displacement = simulation.displacement()
    >>> print(displacement)


.. rst-class:: sphx-glr-script-out

 .. code-block:: none

             results         U
              set_id         3
      node      comp
      4872         X -3.41e-05
                   Y  1.54e-03
                   Z -2.64e-06
      9005         X -5.56e-05
                   Y  1.44e-03
                   Z  5.31e-06
       ...

.. code:: python

    >>> displacement.plot()


.. figure:: ./images/crankshaft_disp.png
    :width: 300pt

.. code:: python

    >>> stress_eqv = simulation.stress_eqv_von_mises_nodal()
    >>> stress_eqv.plot()

.. figure:: ./images/crankshaft_stress.png
    :width: 300pt

To run PyDPF-Post with Ansys 2021 R1 through 2022 R2, use this code to
start the legacy PyDPF-Post tools:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.download_crankshaft())
    >>> stress = solution.stress()
    >>> stress.eqv.plot_contour(show_edges=False)

.. figure:: ./images/crankshaft_stress.png
    :width: 300pt


For comprehensive examples of how you use PyDPF-Post, see :ref:`gallery`.
