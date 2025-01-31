# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from ansys.dpf.core import examples
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.post.faces import Face, FaceListById, FaceListByIndex
from ansys.dpf.post.nodes import NodeListByIndex
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Faces added with ansys-dpf-server 2024.1.pre0.",
)
class TestFaces:
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

    def test_faces(self, fluent_simulation):
        faces = fluent_simulation.mesh.faces
        assert isinstance(faces, FaceListByIndex)
        assert len(faces) == 44242
        ref = "FaceListByIndex([quad4, ..., quad4], __len__=44242)"
        assert repr(faces) == ref
        ref = "[quad4, ..., quad4]"
        assert str(faces) == ref
        assert faces[0] in faces
        for face in faces:
            assert isinstance(face, Face)
        faces_by_id = faces.by_id
        assert isinstance(faces_by_id, FaceListById)
        assert len(faces_by_id) == 44242
        ref = "FaceListById([quad4, ..., quad4], __len__=44242)"
        assert repr(faces_by_id) == ref
        ref = "[quad4, ..., quad4]"
        assert str(faces_by_id) == ref
        assert faces[0] in faces_by_id
        for face in faces_by_id:
            assert isinstance(face, Face)
        assert faces_by_id[faces[0].id].id == faces[0].id
        face = faces[0]
        assert face.node_ids == [11291, 11416, 11455, 11325]
        assert face.id == 1003
        assert face.index == 0
        assert isinstance(face.nodes, NodeListByIndex)
        assert face.num_nodes == 4
        assert face.type == dpf.element_types.Quad4
        ref = "Face(type=element_types.Quad4,index=0,id=1003)"
        assert repr(face) == ref
        ref = """DPF Face 1003
\tIndex:            0
\tNodes:            4
\tType:       element_types.Quad4
"""
        assert str(face) == ref
        assert face.to_node_connectivity == [11290, 11415, 11454, 11324]
