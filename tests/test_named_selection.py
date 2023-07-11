from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf.post import StaticMechanicalSimulation
from ansys.dpf.post.named_selection import NamedSelection


@fixture
def mesh(static_rst):
    simulation = StaticMechanicalSimulation(static_rst)
    return simulation.mesh


def test_named_selections(mesh):
    ns = mesh.named_selections
    assert len(ns) == 1
    assert len(ns.keys()) == 1
    assert all([isinstance(n, str) for n in ns])
    assert ns.keys()[0] == "_FIXEDSU"
    n = ns[ns.keys()[0]]
    assert isinstance(n, NamedSelection)
    assert len(n.ids) == 21


def test_named_selection():
    scoping = dpf.Scoping(ids=[1, 2, 3], location=dpf.locations.nodal)
    ns = NamedSelection(name="test", scoping=scoping)
    assert ns._scoping == scoping
    assert ns.name == "test"
    assert ns.id(0) == 1
    assert ns.index(1) == 0
    assert ns.location == dpf.locations.nodal
