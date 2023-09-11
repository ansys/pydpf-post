from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0
import pytest

from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.tools import streamlines


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
        import pyvista as pv

        pv.OFF_SCREEN = False
        dataframe = simulation.velocity(times=[0.0], zone_ids=[5])
        print(dataframe)
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
