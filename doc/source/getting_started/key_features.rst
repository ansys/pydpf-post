=======================
PyDPF-Post key features
=======================

Computational efficiency
------------------------

PyDPF-Post is based on DPF, whose data framework localizes loading and
postprocessing on the DPF server, enabling rapid postprocessing workflows
because they are written in C and FORTRAN. Because PyDPF-Post presents results
in a Pythonic manner, you can rapidly develop simple or complex postprocessing
scripts.

Easy to use
-----------

The PyDPF-Post API automates the use of chained DPF operators to make postprocessing
easier. By using operators to compute results, you can build your own custom,
low-level scripts to enable fast postprocessing of potentially multi-gigabyte models
with `PyDPF-Core <https://github.com/ansys/pydpf-core>`_.

Ansys solver results files
--------------------------

DPF supports these Ansys solver results files:

- Mechanical APDL (RST, MODE, RFRQ, and RDSP)
- LS-DYNA (D3PLOT and BINOUT)
- Fluent (CAS/DAT.H5 and FLPRJ)
- CFX (CAS/DAT.CFF, FLPRJ, and RES)

For more comprehensive information on file support for Ansys solvers, see the
`main page <https://dpf.docs.pyansys.com/version/stable/index.html>`_
in the PDF-Core documentation.
