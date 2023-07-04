from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Fluid capabilities added with ansys-dpf-server 2024.1.pre0.",
)
class TestSpecies:
    @fixture
    def fluid_simulation(self, fluid_fluent_multi_species):
        ds = dpf.DataSources()
        ds.set_result_file_path(
            fluid_fluent_multi_species["cas"],
            key="cas",
        )
        ds.add_file_path(
            fluid_fluent_multi_species["dat"],
            key="dat",
        )
        return post.FluidSimulation(ds)  # noqa

    def test_species(self):
        species = post.Species(name="test", id=1)
        assert isinstance(species, post.Species)
        assert species.name == "test"
        assert species.id == 1
        assert repr(species) == 'Species<name: "test", id=1>'

    def test_species_list(self, fluid_simulation):
        species_list = post.SpeciesList(simulation=fluid_simulation)
        assert isinstance(species_list, post.SpeciesList)
        assert len(species_list) == 4
        assert isinstance(species_list[0], post.Species)
        ref = (
            '[Species<name: "air", id=0>, Species<name: "ch4", id=1>, '
            'Species<name: "co2", id=2>, Species<name: "o2", id=3>, ]'
        )
        assert repr(species_list) == ref
        ref = "4 species available\n0: air\n1: ch4\n2: co2\n3: o2\n"
        assert str(species_list) == ref
