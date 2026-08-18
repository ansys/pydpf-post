==========
Brief demo
==========

This brief demo shows how to load an Ansys Mechanical APDL (MAPDL) result file to extract
and postprocess results. The code to use depends on which Ansys version you have installed.
For comprehensive examples of how to use PyDPF-Post, see :ref:`gallery`.

2023 R1 and later
-----------------

If Ansys 2023 R1 or later is installed, a DPF server automatically starts
once you start using PyDPF-Post. Use this code to load an MAPDL result file
to extract and postprocess results:

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


.. image:: ./../images/crankshaft_disp.png
    :align: center
    :width: 300pt


.. code:: python

    >>> stress_eqv = simulation.stress_eqv_von_mises_nodal()
    >>> stress_eqv.plot()

.. image:: ./../images/crankshaft_stress.png
    :align: center
    :figwidth: 300pt


2021 R1 through 2022 R2
-----------------------

If an Ansys release of 2021 R1 through 2022 R2 is installed, use this code to
start the legacy PyDPF-Post tools and then load an MAPDL result file
to extract and postprocess results:

.. code:: python

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.download_crankshaft())
    >>> stress = solution.stress()
    >>> stress.eqv.plot_contour(show_edges=False)

.. image:: ./../images/crankshaft_stress.png
    :align: center
    :width: 300pt

