"""PyDPF-Post example files.

Examples
--------
This module module exposes PyDPF-Core functionalities.

See `here <https://dpf.docs.pyansys.com/dev/api/ansys.dpf.core.examples.html>`_
for a description of the PyDPF-Core ``example`` API.

"""
import os

# alias files used by the core
from ansys.dpf.core.examples import *  # noqa: F403

# must copy files over when using docker
# running_docker = os.environ.get('DPF_DOCKER', False)
