from ansys.dpf.core import examples
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
class TestPhase:
    @fixture
    def fluid_simulation(self):
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

    def test_phase(self):
        phase = post.Phase(name="phase_test", id=1)
        assert isinstance(phase, post.Phase)
        assert phase.name == "phase_test"
        assert phase.id == 1
        assert repr(phase) == 'Phase<name: "phase_test", id=1>'

    def test_phases(self, fluid_simulation):
        phases = post.Phases(simulation=fluid_simulation)
        assert isinstance(phases, post.Phases)
        assert len(phases) == 1
        assert isinstance(phases[1], post.Phase)
        assert repr(phases) == "[Phase<name: 'phase-1', id=1>, ]"
        print(phases)
        assert str(phases) == "1 phases available\n{Phase<name: 'phase-1', id=1>,}"
        with pytest.raises(ValueError):
            _ = phases["toto"]
        with pytest.raises(ValueError):
            _ = phases[32]
        with pytest.raises(ValueError):
            _ = phases[24.6]
        for phase in phases:
            assert isinstance(phase, post.Phase)
        for phase in phases:
            assert isinstance(phase, post.Phase)
        _ = list(phases)
