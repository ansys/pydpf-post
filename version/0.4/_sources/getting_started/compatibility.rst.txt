.. _compatibility:

=============
Compatibility
=============

PyDPF-Post supports Windows 10 and CentOS 7 and later. For
more information, see `Ansys Platform Support <https://www.ansys.com/solutions/solutions-by-role/it-professionals/platform-support>`_.

Other platforms may be supported by using PyDPF-Post in a
containerized ecosystem, such as `Docker <https://www.docker.com/>`_ or
`Kubernetes <https://kubernetes.io/>`_.

Due to potential crashes when using PyDPF-Post 0.2.2 or earlier with PyDPF-Core 0.5.2
or later, you should refer to the following table for compatibility between these two
libraries.

.. list-table:: PyDPF compatibility
   :widths: 20 20 20
   :header-rows: 1

   * - DPF server version
     - ansys.dpf.core python module version
     - ansys.dpf.post python module version
   * - 6.1 (DPF Server 2023.2.pre1)
     - 0.8.0 or later
     - 0.4.0 or later
   * - 6.0 (DPF Server 2023.2.pre0)
     - 0.7.0 or later
     - 0.3.0 or later
   * - 5.0 (Ansys 2023 R1)
     - 0.6.0 or later
     - 0.3.0 or later
   * - 4.0 (Ansys 2022 R2)
     - 0.5.0 or later
     - 0.3.0 or later
   * - 3.0 (Ansys 2022 R1)
     - 0.4.0 or later
     - 0.1.0 or later
   * - 2.0 (Ansys 2021 R2)
     - 0.3.0 or later
     - 0.1.0 or later
   * - 1.0 (Ansys 2021 R1)
     - 0.2.*
     - 0.1.0 or later
