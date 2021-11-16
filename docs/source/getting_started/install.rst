.. _installation:

************
Installation
************
Once you've installed Ansys 2021R1 or newer, you can install DPF with:

.. code::

   pip install ansys-dpf-post


This will install the latest version of ``ansys-dpf-post`` and all the
necessary dependencies.

If you are unable to install the module on the host machine due to
network isolation, download the latest release wheel at `pydpf-post
GitHub <https://github.com/pyansys/pydpf-post>`_ or from PyPi at
`pydpf-post PyPi <https://pypi.org/project/ansys-dpf-post/>`_


Editable Install (Development Mode)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you wish to edit and potentially contribute to the DPF-Post python
module, clone the repository and install it using pip with the ``-e``
development flag.

.. code::

    git clone https://github.com/pyansys/DPF-Post
    cd DPF-Post
    pip install -e .
