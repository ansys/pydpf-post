.. _user_guide_accessing_file_metadata:

**********************
Browse result metadata
**********************

In addition to the :class:`DpfSolution <ansys.dpf.post.dpf_solution.DpfSolution>`
object being the entry point for browsing the contents of a result file, it provides
important metadata, such as the analysis type and the available results.

Here is how you browse result metadata:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Browse result metadata

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


DPF-Post supports four different analysis types:

* Static analysis
* Modal analysis
* Harmonic analysis
* Transient analysis

Mesh
----
From the ``Solution`` object, you can also access the mesh:

.. code:: python

    Instantiate the solution object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.multishells_rst)

    Access the mesh

    >>> mesh = solution.mesh
    >>> print(mesh)
    Meshed Region
    	7079 nodes
    	4220 elements
    	Unit: m

