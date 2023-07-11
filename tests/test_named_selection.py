import pytest
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
    with pytest.raises(KeyError, match="could not be found"):
        _ = ns["test"]


def test_named_selection():
    scoping = dpf.Scoping(ids=[1, 2, 3], location=dpf.locations.nodal)
    ns = NamedSelection(name="test", scoping=scoping)
    assert ns._scoping == scoping
    assert ns.name == "test"
    assert ns.id(0) == 1
    assert ns.index(1) == 0
    assert ns.location == dpf.locations.nodal
    assert ns.size == 3
    ref = """NamedSelection 'test'
 with DPF  Scoping: 
  with Nodal location and 3 entities
"""  # noqa
    assert repr(ns) == ref

    assert ns.deep_copy() == ns

    assert ns.as_local_scoping() == ns

    with pytest.raises(ValueError, match="No list of entity IDs given."):
        _ = NamedSelection(name="test")
    with pytest.raises(ValueError, match="NamedSelection accepts only"):
        _ = NamedSelection(name="test", node_ids=[0, 1], element_ids=[0, 1])

    node_ns = NamedSelection(name="test", node_ids=[0, 1])
    assert node_ns._scoping.location == dpf.locations.nodal
    element_ns = NamedSelection(name="test", element_ids=[0, 1])
    assert element_ns._scoping.location == dpf.locations.elemental
    face_ns = NamedSelection(name="test", face_ids=[0, 1])
    assert face_ns._scoping.location == dpf.locations.faces
    cell_ns = NamedSelection(name="test", cell_ids=[0, 1])
    assert cell_ns._scoping.location == dpf.locations.elemental
