.. _user_guide_accessing_file_metadata:

**********************
Browse result metadata
**********************

In addition to the :class:`Simulation <ansys.dpf.post.simulation.Simulation>`
object being the entry point for browsing the contents of a result file, it provides
important metadata, such as the analysis type and the available results.

Here is how you browse result metadata:

.. code:: python

    Instantiate the simulation object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> simulation = post.load_simulation(examples.multishells_rst)

    Browse result metadata

    >>> print(simulation)
    Static Mechanical Simulation.


    Data Sources
    ------------------------------
    DPF  DataSources:
      Result files:
         result key: rst and path: d:\ansysdev\sandbox\pydpf-core\src\ansys\dpf\core\examples\model_with_ns.rst
      Secondary files:


    DPF Model
    ------------------------------
    Static analysis
    Unit system: MKS: m, kg, N, s, V, A, degC
    Physics Type: Mechanical
    Available results:
         -  displacement: Nodal Displacement
         -  reaction_force: Nodal Force
         -  stress: ElementalNodal Stress
         -  elemental_volume: Elemental Volume
         -  stiffness_matrix_energy: Elemental Energy-stiffness matrix
         -  artificial_hourglass_energy: Elemental Hourglass Energy
         -  thermal_dissipation_energy: Elemental thermal dissipation energy
         -  kinetic_energy: Elemental Kinetic Energy
         -  co_energy: Elemental co-energy
         -  incremental_energy: Elemental incremental energy
         -  elastic_strain: ElementalNodal Strain
         -  structural_temperature: ElementalNodal Temperature
    ------------------------------
    DPF  Meshed Region:
      7079 nodes
      4220 elements
      Unit: m
      With solid (3D) elements, shell (2D) elements, shell (3D) elements
    ------------------------------
    DPF  Time/Freq Support:
      Number of sets: 1
    Cumulative     Time (s)       LoadStep       Substep
    1              1.000000       1              1


PyDPF-Post supports four different analysis types for mechanical results:

* Static analysis
* Modal analysis
* Harmonic analysis
* Transient analysis

The legacy Solution object also supports thermal and electrical results.

Mesh
----
From the ``Simulation`` object, you can also access the mesh:

.. code:: python

    Instantiate the simulation object

    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> simulation = post.load_simulation(examples.multishells_rst)

    Access the mesh

    >>> mesh = simulation.mesh
    >>> print(mesh)
    DPF  Mesh:
      7079 nodes
      4220 elements
      Unit: m
      With solid (3D) elements, shell (2D) elements, shell (3D) elements

