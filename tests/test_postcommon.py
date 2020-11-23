from ansys import dpf
from ansys.dpf import post


if not dpf.core.has_local_server():
    dpf.core.start_local_server()
    
    
def test_available_keywords():
    txt = post.available_keywords().__str__()
    assert "el_shape" in txt
    assert "element_scoping" in txt
    assert "grouping" in txt
    assert "location" in txt
    assert "named_selection" in txt
    assert "node_scoping" in txt
    assert "phase" in txt