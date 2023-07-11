import pytest

from ansys.dpf import core as dpf
from ansys.dpf.post.elements import Element, ElementType


def test_element_type():
    with pytest.raises(TypeError, match="Given argument"):
        _ = ElementType("test")

    el_type = ElementType(
        arg=dpf.ElementDescriptor(
            enum_id=dpf.element_types.Beam4, description="test", name="test"
        )
    )
    ref = "test"
    assert el_type.description == ref

    el_type = ElementType(arg=1)
    ref = """Element Type
------------
Enum id (dpf.element_types): element_types.Hex20
Element description: Quadratic 20-nodes Hexa
Element name (short): hex20
Element shape: solid
Number of corner nodes: 8
Number of mid-side nodes: 12
Total number of nodes: 20
Quadratic element: True"""
    assert str(el_type) == ref
    assert repr(el_type) == ref
    ref = "Quadratic 20-nodes Hexa"
    assert el_type.description == ref
    assert el_type.enum_id == dpf.element_types.Hex20
    assert el_type.name == "hex20"
    assert el_type.shape == "solid"
    assert el_type.num_corner_nodes == 8
    assert el_type.num_mid_nodes == 12
    assert el_type.num_nodes == 20
    assert el_type.is_quadratic
    assert el_type.is_solid
    assert not el_type.is_shell
    assert not el_type.is_beam


def test_element():
    element = Element()
