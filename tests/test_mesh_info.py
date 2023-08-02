from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.post import examples


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="FluidMeshInfo added with ansys-dpf-server 2024.1.pre0.",
)
class TestFluidMeshInfo:
    @fixture
    def fluent_simulation(self):
        fluid_example_files = examples.download_fluent_axial_comp()
        ds = dpf.DataSources()
        ds.set_result_file_path(
            fluid_example_files["cas"][0],
            key="cas",
        )
        ds.add_file_path(
            fluid_example_files["dat"][0],
            key="dat",
        )
        return post.FluidSimulation(ds)  # noqa

    def test_fluid_mesh_info_print(self, fluent_simulation):
        # print(fluent_simulation.mesh_info)
        ref = (
            "Fluid mesh metadata\n"
            "-------------------\n"
            "Number of nodes: 16660\n"
            "Number of faces: 45391\n"
            "Number of cells: 13856\n"
            "Cell zones:\n"
            "\t{13: 'fluid-rotor', 28: 'fluid-stator'}\n"
            "Face zones:\n"
            "\t{2: 'default-interior:0', 3: 'rotor-hub', 4: 'rotor-shroud', 5: 'rotor-inlet', "
            "6: 'rotor-interface', 7: 'rotor-blade-1', 8: 'rotor-blade-2', "
            "9: 'rotor-per-1-shadow', 10: 'rotor-per-1', 11: 'rotor-per-2-shadow', "
            "12: 'rotor-per-2', 15: 'default-interior', 16: 'stator-hub', "
            "17: 'stator-shroud', 18: 'stator-interface', 19: 'stator-outlet', "
            "20: 'stator-blade-1', 21: 'stator-blade-2', 22: 'stator-blade-3', "
            "23: 'stator-blade-4', 24: 'stator-per-2', 25: 'stator-per-2-shadow', "
            "26: 'stator-per-1', 27: 'stator-per-1-shadow'}\n"
            "Cell to face zones:\n"
            "\t{13: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "
            "28: [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]}\n"
        )
        assert str(fluent_simulation.mesh_info) == ref
