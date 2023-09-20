import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0


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

    def test_species_dict(self, fluid_simulation):
        species_dict = post.SpeciesDict(simulation=fluid_simulation)
        assert isinstance(species_dict, post.SpeciesDict)
        assert len(species_dict) == 4
        assert isinstance(species_dict[1], post.Species)
        ref = (
            '{Species<name: "air", id=4>, Species<name: "ch4", id=3>, '
            'Species<name: "co2", id=1>, Species<name: "o2", id=2>, }'
        )
        assert repr(species_dict) == ref
        ref = (
            "4 species available\n"
            '{Species<name: "air", id=4>, Species<name: "ch4", id=3>, '
            'Species<name: "co2", id=1>, Species<name: "o2", id=2>, }'
        )
        assert str(species_dict) == ref
        with pytest.raises(ValueError):
            _ = species_dict["toto"]
        with pytest.raises(ValueError):
            _ = species_dict[32]
        with pytest.raises(ValueError):
            _ = species_dict[24.6]
        for species in species_dict:
            assert isinstance(species, post.Species)
        for species in species_dict:
            assert isinstance(species, post.Species)
        _ = list(species_dict)
