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

import pytest

from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.helpers import streamlines
from conftest import (
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0,
)


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
class TestStreamlines:
    @pytest.fixture
    def simulation(self) -> post.FluidSimulation:
        files_cfx = examples.download_cfx_heating_coil()
        return post.FluidSimulation(cas=files_cfx["cas"], dat=files_cfx["dat"])  # noqa

    def test_plot_streamlines(self, simulation):
        if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0:
            dataframe = simulation.velocity(times=[0.0], zone_ids=[1])
        else:
            dataframe = simulation.velocity(times=[0.0], zone_ids=[5])
        sources = [
            {"radius": 0.25, "center": (0.75, 0.0, 0.0), "n_points": 20},
            {
                "radius": 0.25,
                "center": (0.0, 0.75, 0.0),
                "n_points": 5,
                "max_time": 10.0,
            },
            {"radius": 0.25, "center": (-0.75, 0.0, 0.0), "max_time": 2.0},
            {"radius": 0.25, "center": (0.0, -0.75, 0.0)},
        ]
        streamlines.plot_streamlines(
            dataframe=dataframe,
            sources=sources,
            streamline_thickness=0.007,
            plot_mesh=True,
            mesh_opacity=0.2,
            plot_contour=True,
            contour_opacity=0.3,
            title="Streamlines with multiple sources",
        )
        dataframe = simulation.velocity()
        with pytest.raises(ValueError, match="The set_id requested is not available"):
            streamlines.plot_streamlines(
                dataframe=dataframe,
                sources=sources,
                set_id=2,
                streamline_thickness=0.007,
                plot_mesh=True,
                mesh_opacity=0.2,
                plot_contour=True,
                contour_opacity=0.3,
                title="Streamlines with multiple sources",
            )
        streamlines.plot_streamlines(
            dataframe=dataframe,
            sources=sources,
            set_id=1,
            streamline_thickness=0.007,
            plot_mesh=True,
            mesh_opacity=0.2,
            plot_contour=True,
            contour_opacity=0.3,
            title="Streamlines with multiple sources",
        )
