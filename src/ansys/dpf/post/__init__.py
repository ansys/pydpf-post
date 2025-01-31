# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Post-processing module.

Using the Data Processing Framework.

Allow to create a result object, then use it to get desired results.

Examples
--------
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> solution = post.load_solution(examples.static_rst)
>>> disp = solution.displacement()

"""

import ansys.dpf.core as core
from ansys.dpf.core import (  # noqa: F401
    AvailableServerContexts,
    set_default_server_context,
)
from ansys.dpf.core.common import locations, shell_layers  # noqa: F401

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

from ansys.dpf.post import mesh, selection, tools  # noqa: F401
from ansys.dpf.post.common import Grouping as grouping  # noqa: F401
from ansys.dpf.post.dataframe import DataFrame  # noqa: F401
from ansys.dpf.post.dpf_path import create_path_on_coordinates  # noqa: F401
from ansys.dpf.post.fluid_simulation import FluidSimulation  # noqa: F401
from ansys.dpf.post.harmonic_mechanical_simulation import (  # noqa: F401
    HarmonicMechanicalSimulation,
)
from ansys.dpf.post.mesh import Mesh  # noqa: F401
from ansys.dpf.post.meshes import Meshes  # noqa: F401
from ansys.dpf.post.misc import Report
from ansys.dpf.post.modal_mechanical_simulation import (  # noqa: F401
    ModalMechanicalSimulation,
)
from ansys.dpf.post.phase import Phase, PhasesDict  # noqa: F401
from ansys.dpf.post.post_utility import (  # noqa: F401
    load_simulation,
    load_solution,
    print_available_keywords,
)
from ansys.dpf.post.species import Species, SpeciesDict  # noqa: F401
from ansys.dpf.post.static_mechanical_simulation import (  # noqa: F401
    StaticMechanicalSimulation,
)
from ansys.dpf.post.transient_mechanical_simulation import (  # noqa: F401
    TransientMechanicalSimulation,
)

# this must be after some ansys.dpf.post import
__version__ = importlib_metadata.version("ansys-dpf-post")
