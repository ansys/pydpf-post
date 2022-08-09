===============
Getting Started
===============
To use the PyDPF-Post API, you need to have the PyDPF-Core module installed.
You also need to use a native DPF server and thus need a local installation of
Ansys 2021R1 or newer.  Visit https://www.ansys.com/ for more information on
getting a licensed copy of Ansys.


.. toctree::
   :hidden:
   :maxdepth: 2

   install


Compatibility
~~~~~~~~~~~~~
PyDPF-Post supports Windows 10 and CentOS 7 and newer. For
more details, see `Ansys Platform Support <https://www.ansys.com/solutions/solutions-by-role/it-professionals/platform-support>`_.

Other platforms may be supported by using PyDPF-Post within a
containerization ecosystem such as Docker or Kubernetes.

Due to potential crashes when using PyDPF-Post <=0.2.2 with PyDPF-Core >=0.5.2,
it is advised to follow the table below for compatibility between PyDPF modules.

.. list-table:: PyDPF Compatibility
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Ansys installation
     - ansys.dpf.core python module version
     - ansys.dpf.post python module version
   * - Ansys 2022R2
     - >=0.5.0
     - >=0.3.0
   * - Ansys 2022R1
     - >=0.4.0
     - >=0.1.0
   * - Ansys 2021R2
     - >=0.3.0
     - >=0.1.0
   * - Ansys 2021R1
     - 0.2.*
     - >=0.1.0