from ansys import dpf
from ansys.dpf import post


def test_available_keywords():
    txt = post.common._AvailableKeywords().__str__()
    assert "element_scoping" in txt
    assert "grouping" in txt
    assert "location" in txt
    assert "named_selection" in txt
    assert "node_scoping" in txt
