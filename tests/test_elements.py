from ansys.dpf.core.nodes import Node
import pytest

from ansys.dpf import core as dpf
from ansys.dpf.post import examples
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
    model = dpf.Model(examples.find_static_rst())
    core_elements = model.metadata.meshed_region.elements
    element = Element(element=core_elements[0])
    ref = [1, 26, 14, 12, 2, 27, 15, 13, 33, 64, 59, 30, 37, 65, 61, 34, 28, 81, 63, 58]
    assert element.node_ids == ref
    assert element.id == 5
    assert element.index == 0
    assert isinstance(element.nodes[0], Node)
    assert element.num_nodes == 20
    assert isinstance(element.type_info, ElementType)
    assert element.type == dpf.element_types.Hex20
    assert element.shape == "solid"
    ref = [0, 25, 13, 11, 1, 26, 14, 12, 32, 63, 58, 29, 36, 64, 60, 33, 27, 80, 62, 57]
    assert element.to_node_connectivity == ref
    ref = """DPF Element 5
\tIndex:            0
\tNodes:           20
\tType:         Hex20
\tShape:        Solid
"""
    assert str(element) == ref
    ref = "Element(type=element_types.Hex20,index=0,id=5,shape=solid)"
    assert repr(element) == ref
