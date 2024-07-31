import numpy as np
import pytest
from pytest import fixture

from ansys.dpf import core as dpf
from ansys.dpf import post
from ansys.dpf.post import examples
from ansys.dpf.post.selection import SpatialSelection
from conftest import SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0


def test_spatial_selection_select_nodes(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection._selection.progress_bar = False
    selection.select_nodes([1, 2, 3])
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.locations.nodal
    assert np.allclose(scoping.ids, [1, 2, 3])


def test_spatial_selection_select_elements(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection._selection.progress_bar = False
    selection.select_elements([1, 2, 3, 4])
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.locations.elemental
    assert np.allclose(scoping.ids, [1, 2, 3, 4])
    ids = selection.apply_to(simulation)
    assert np.allclose(ids, [1, 2, 3, 4])


def test_spatial_selection_select_named_selection(allkindofcomplexity):
    simulation = post.load_simulation(allkindofcomplexity)
    selection = SpatialSelection()
    selection._selection.progress_bar = False
    selection.select_named_selection(
        simulation.mesh.named_selections.keys()[0],
        location=post.selection.locations.nodal,
    )
    scoping = selection._evaluate_on(simulation)
    assert scoping.location == post.locations.nodal
    assert scoping.ids.size == 12970
    assert 1857 in scoping.ids
    assert 14826 in scoping.ids
    ids = selection.apply_to(simulation)
    assert len(ids) == 12970
    assert 1857 in ids
    assert 14826 in ids


@pytest.mark.skipif(
    not SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0,
    reason="Faces added with ansys-dpf-server 2024.1.pre0.",
)
class TestSpatialSelectionFaces:
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

    def test_spatial_selection_select_faces(self, fluent_simulation):
        selection = SpatialSelection()
        selection._selection.progress_bar = False
        selection.select_faces(fluent_simulation.mesh.face_ids)
        scoping = selection._evaluate_on(fluent_simulation)
        assert scoping.location == post.locations.faces
        assert np.allclose(scoping.ids, fluent_simulation.mesh.face_ids)

    def test_spatial_selection_select_nodes_of_faces(self, fluent_simulation):
        selection = SpatialSelection()
        selection._selection.progress_bar = False
        face_0 = fluent_simulation.mesh.faces[0]
        selection.select_nodes_of_faces(
            faces=[face_0.id],
            mesh=fluent_simulation.mesh,
        )
        scoping = selection._evaluate_on(fluent_simulation)
        assert scoping.location == post.locations.nodal
        assert np.allclose(scoping.ids, face_0.node_ids)

    def test_spatial_selection_select_faces_of_elements(self, fluent_simulation):
        selection = SpatialSelection()
        selection._selection.progress_bar = False
        elem_0 = fluent_simulation.mesh.elements[0]
        selection.select_faces_of_elements(
            elements=[elem_0.id],
            mesh=fluent_simulation.mesh,
        )
        scoping = selection._evaluate_on(fluent_simulation)
        assert scoping.location == post.locations.faces
        assert np.allclose(scoping.ids, [11479, 11500, -1, 11502, 11503])


    def test_skin_extraction_with_scoping(self):

        # Todo: Upload to examples repo
        rst_path = r"D:\ANSYSDev\remote_post\models\beam_2_mat.rst"

        simulation = post.StaticMechanicalSimulation(rst_path)

        # get an elemental named selection
        scoping = simulation._model.metadata.meshed_region.named_selection(
            'X_MIN_SOLID_ELE'
        )

        seqv_ns_skin = simulation.stress_eqv_von_mises_elemental(skin=scoping.ids)

        field = seqv_ns_skin._fc[0]

        assert field.scoping.size == 128
        assert field.max().data[0] == pytest.approx(287.2531)
        assert field.min().data[0] == pytest.approx(213.1394)

#
#
# def test_spatial_selection_intersect(allkindofcomplexity):
#     solution = post.load_solution(allkindofcomplexity, legacy=False)
#     selection1 = SpatialSelection()
#     selection1.select_nodes([1, 2, 3])
#     _ = selection1._evaluate_on(solution)
#
#     selection = SpatialSelection()
#     selection.select_nodes([1, 2, 3, 4])
#     selection.intersect(selection1)
#     scoping = selection._evaluate_on(solution)
#     assert scoping.location == post.selection.locations.nodal
#     assert np.allclose(scoping.ids, [1, 2, 3])
#     ids = selection.apply_to(solution)
#     assert np.allclose(ids, [1, 2, 3])
