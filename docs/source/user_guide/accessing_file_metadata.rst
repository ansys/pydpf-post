.. _user_guide_accessing_file_metadata:

************************
Browsing Result Metadata
************************

The ``Solution`` object is the entry point for browsing the contents
of a result file (see :ref:`user_guide_post_processing`). It provides
access to important metadata such as the analysis type and listings of
available results.

This information is accessible as follows.

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Access the result information

    >>> result_info = solution.get_result_info()
    >>> print(result_info)
    Static analysis
    Unit system: Metric (m, kg, N, s, V, A)
    Physics Type: Mecanic
    Available results:
         -  displacement
         -  force
         -  stress
         -  volume
         -  energy_stiffness_matrix
         -  hourglass_energy
         -  thermal_dissipation_energy
         -  kinetic_energy
         -  co_energy
         -  incremental_energy
         -  strain
         -  temperature


Four different analysis types are presently supported in DPF-Post:

* static analysis
* modal analysis
* harmonic analysis
* transient analysis

Mesh
----
The mesh can also be accessed from the ``Solution`` object.

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Get the mesh

    >>> mesh = solution.mesh
    >>> print(mesh)
    Meshed Region
    	7079 nodes
    	4220 elements
    	Unit: m
