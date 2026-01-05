# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
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

"""
.. _ref_cyclic_mesh_skin_example:

Reduce cyclic model size by using the mesh skin for result and mesh extraction
==============================================================================
This example shows postprocessing on a mesh skin for a cyclic model analysis.
The skin mesh is rebuilt with new surface elements connecting the nodes on the external skin
of the solid mesh. These surface elements types are chosen with respect to the solid elements
facets having all their nodes on the skin.

This feature, available for all types of mechanical simulation supporting cyclic
or cyclic multi-stage models, allows you to reduce the size of both the mesh
and extracted data to improve processing performance. Because larger stresses and
strains are usually located on the skin of a model, computing results on the skin gives
equivalent maximum values in most cases.

Postprocessing of elemental or elemental nodal results requires an element solid-to-skin mapping
to get from a solid element result to a facet result. Because the connectivity of the new surface
elements built on the skin are different from the connectivity of the solid elements, small
differences can be found after result averaging.

To plot cyclic expanded results, the skin mesh is expanded.
"""

###############################################################################
# Perform required imports
# ------------------------
# Perform required imports. This example uses a supplied file that you can
# get using the ``examples`` module.

from ansys.dpf import post
from ansys.dpf.post import examples

###############################################################################
# Get ``Simulation`` object
# -------------------------
# Get the ``Simulation`` object that allows access to the result. The ``Simulation``
# object must be instantiated with the path for the result file. For example,
# ``"C:/Users/user/my_result.rst"`` on Windows or ``"/home/user/my_result.rst"``
# on Linux.

example_path = examples.download_modal_cyclic()
simulation = post.ModalMechanicalSimulation(example_path)

# print the simulation to get an overview of what's available
print(simulation)


###############################################################################
# Extract displacement data
# -------------------------
# Extract displacement data on the skin.

displacement_skin = simulation.displacement(skin=True)
displacement_skin.plot()


###############################################################################
# Extract stress and strain data
# ------------------------------
# Extract stress and elastic strain data over the entire mesh and on the skin.
# Averaging and invariants computation are done through a solid-to-skin connectivity mapping.

elemental_stress_skin = simulation.stress_principal_elemental(components=[1], skin=True)
elemental_stress_skin.plot()

elastic_strain_eqv_skin = simulation.elastic_strain_eqv_von_mises_nodal(skin=True)
elastic_strain_eqv_skin.plot()

###############################################################################
# Get stress results on skin of first sector with a cyclic phase
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stress_eqv_cyc_phase = simulation.stress_eqv_von_mises_nodal(
    set_ids=[5],
    expand_cyclic=[1],
    phase_angle_cyclic=45.0,
    skin=True,
)
stress_eqv_cyc_phase.plot()
