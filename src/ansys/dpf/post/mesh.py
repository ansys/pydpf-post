"""Module containing the ``Mesh`` class.

Mesh
----

"""
from __future__ import annotations

from typing import List

import ansys.dpf.core as dpf
from ansys.dpf.core.nodes import Node
from ansys.dpf.core.property_fields_container import (
    _MockPropertyFieldsContainer as PropertyFieldsContainer,
)

import ansys.dpf.post as post
from ansys.dpf.post import index, locations
from ansys.dpf.post.connectivity import ConnectivityListByIndex, ReturnMode
from ansys.dpf.post.elements import Element, ElementListByIndex
from ansys.dpf.post.faces import FaceListByIndex
from ansys.dpf.post.named_selection import NamedSelections
from ansys.dpf.post.nodes import NodeListByIndex


class Mesh:
    """Exposes the complete mesh of the simulation."""

    def __init__(self, meshed_region: dpf.MeshedRegion):
        """Initialize this class."""
        if meshed_region is None:
            raise ValueError("Tried to instantiate an empty Mesh.")
        self._meshed_region = meshed_region

    def __str__(self):
        """String representation of this class."""
        return str(self._meshed_region).replace("Meshed Region", "Mesh")

    @property
    def named_selections(self) -> NamedSelections:
        """Returns a dictionary of available named selections for this mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.named_selections) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        NamedSelections dictionary with 1 named selections:
            - '_FIXEDSU'
        """
        return NamedSelections(self)

    @property
    def node_ids(self) -> List[int]:
        """Returns the list of node IDs in the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.node_ids) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
         25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
         49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72
         73 74 75 76 77 78 79 80 81]
        """
        return self._meshed_region.nodes.scoping.ids

    @property
    def num_nodes(self) -> int:
        """Returns the number of nodes in the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.num_nodes) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        81
        """
        return self._meshed_region.nodes.n_nodes

    @property
    def element_ids(self) -> List[int]:
        """Returns the list of element IDs in the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.element_ids) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [5 6 1 2 7 8 3 4]
        """
        return self._meshed_region.elements.scoping.ids

    @property
    def num_elements(self) -> int:
        """Returns the number of elements in the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.num_elements) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        8
        """
        return self._meshed_region.elements.n_elements

    @property
    def elements(self) -> ElementListByIndex:
        """Returns a list of elements indexed by ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.elements) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [hex20, ..., hex20]
        """
        return ElementListByIndex(self._meshed_region.elements)

    def get_element_by_id(
        self, id: int  # pylint: disable=redefined-builtin
    ) -> Element:
        """Returns an element in the mesh from a given ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.get_element_by_id(1)) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        DPF Element 1
            Index:            2
            Nodes:           20
            Type:         Hex20
            Shape:        Solid
        """
        return self.elements.by_id[id]

    @property
    def num_faces(self) -> int:
        """Returns the number of faces in the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> files = examples.download_fluent_axial_comp()
        >>> simulation = post.FluidSimulation(cas=files['cas'][0], dat=files['dat'][0])
        >>> print(simulation.mesh.num_faces) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        44242
        """
        return self._meshed_region.faces.n_faces

    @property
    def face_ids(self) -> List[int]:
        """Returns the list of face IDs in the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> files = examples.download_fluent_axial_comp()
        >>> simulation = post.FluidSimulation(cas=files['cas'][0], dat=files['dat'][0])
        >>> print(simulation.mesh.face_ids) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [ 1003  1004  1005 ... 45165 45166 45167]
        """
        return self._meshed_region.faces.scoping.ids

    @property
    def faces(self) -> FaceListByIndex:
        """Returns a list of faces indexed by ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> files = examples.download_fluent_axial_comp()
        >>> simulation = post.FluidSimulation(cas=files['cas'][0], dat=files['dat'][0])
        >>> print(simulation.mesh.faces) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [quad4, ..., quad4]
        """
        return FaceListByIndex(self._meshed_region.faces)

    def get_face_by_id(self, id: int) -> Face:  # pylint: disable=redefined-builtin
        """Returns a face in the mesh from a given ID."""
        return self.faces.by_id[id]

    @property
    def nodes(self) -> NodeListByIndex:
        """Returns a list of nodes indexed by ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.nodes) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [Node(id=1, coordinates=[0.015, 0.045, 0.015]), ..., Node(id=81, coordinates=[0.03, 0.045, 0.0075])]
        """  # noqa
        return NodeListByIndex(self._meshed_region.nodes)

    def get_node_by_id(self, id: int) -> Node:  # pylint: disable=redefined-builtin
        """Returns a node of the mesh from a given ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.get_node_by_id(1)) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        Node(id=1, coordinates=[0.015, 0.045, 0.015])
        """
        return self.nodes.by_id[id]

    @property
    def element_types(self) -> post.DataFrame:
        """Returns a DataFrame containing element types ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.element_types) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
             results elem_type_id
         element_ids
                   5            1
                   6            1
                   1            1
                   2            1
                   7            1
                   8            1
                 ...          ...
        """
        label = "elem_type_id"
        fields_container = PropertyFieldsContainer()
        field = self._meshed_region.elements.element_types_field
        fields_container.add_field(label_space={}, field=field)

        return post.DataFrame(
            data=fields_container,
            index=index.MultiIndex(
                indexes=[
                    index.MeshIndex(
                        location=locations.elemental,
                        scoping=self._meshed_region.elements.scoping,
                        fc=fields_container,
                    )
                ]
            ),
            columns=index.MultiIndex(indexes=[index.ResultsIndex(values=[label])]),
        )

    @property
    def materials(self) -> post.DataFrame:
        """Returns a DataFrame containing element materials ID.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.materials) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
             results material_id
         element_ids
                   5           1
                   6           1
                   1           1
                   2           1
                   7           1
                   8           1
                 ...         ...
        """
        label = "material_id"
        fields_container = PropertyFieldsContainer()
        field = self._meshed_region.elements.materials_field
        fields_container.add_field(label_space={}, field=field)

        return post.DataFrame(
            data=fields_container,
            index=index.MultiIndex(
                indexes=[
                    index.MeshIndex(
                        location=locations.elemental,
                        scoping=self._meshed_region.elements.scoping,
                        fc=fields_container,
                    )
                ]
            ),
            columns=index.MultiIndex(indexes=[index.ResultsIndex(values=[label])]),
        )

    @property
    def element_to_node_ids_connectivity(self) -> ConnectivityListByIndex:
        """Returns a connectivity map between element index and node IDs.

        To get the connectivity map by element ID, use the 'by_id' property of the object returned.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> element_to_node_ids_connectivity = simulation.mesh.element_to_node_ids_connectivity
        >>> print(element_to_node_ids_connectivity[0]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [1, 26, 14, 12, 2, 27, 15, 13, 33, 64, 59, 30, 37, 65, 61, 34, 28, 81, 63, 58]
        >>> element_to_node_ids_connectivity_by_id = element_to_node_ids_connectivity.by_id
        >>> print(element_to_node_ids_connectivity[1]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [1, 12, 14, 26, 3, 10, 9, 4, 30, 59, 64, 33, 41, 53, 43, 38, 29, 56, 54, 44]
        """
        conn_field = self._meshed_region.elements.connectivities_field
        nodes_scoping = self._meshed_region.nodes.scoping
        return ConnectivityListByIndex(
            field=conn_field, mode=ReturnMode.IDS, scoping=nodes_scoping
        )

    @property
    def node_to_element_ids_connectivity(self) -> ConnectivityListByIndex:
        """Returns a connectivity map between node index and element IDs.

        To get the connectivity map by node ID, use the 'by_id' property of the object returned.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> node_to_element_ids_connectivity = simulation.mesh.node_to_element_ids_connectivity
        >>> print(node_to_element_ids_connectivity[0]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [5, 6, 1, 2, 7, 8, 3, 4]
        >>> node_to_element_ids_connectivity_by_id = node_to_element_ids_connectivity.by_id
        >>> print(node_to_element_ids_connectivity_by_id[1]) # doctest: +NORMALIZE_WHITESPACE
        [5, 6, 1, 2, 7, 8, 3, 4]
        """
        conn_field = self._meshed_region.nodes.nodal_connectivity_field
        elems_scoping = self._meshed_region.elements.scoping
        return ConnectivityListByIndex(
            field=conn_field, mode=ReturnMode.IDS, scoping=elems_scoping
        )

    @property
    def element_to_node_connectivity(self) -> ConnectivityListByIndex:
        """Returns a connectivity map between element index and node indexes.

        To get the connectivity map by element ID, use the 'by_id' property of the object returned.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> element_to_node_connectivity = simulation.mesh.element_to_node_connectivity
        >>> print(element_to_node_connectivity[0]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [0, 25, 13, 11, 1, 26, 14, 12, 32, 63, 58, 29, 36, 64, 60, 33, 27, 80, 62, 57]
        >>> element_to_node_connectivity_by_id = element_to_node_connectivity.by_id
        >>> print(element_to_node_connectivity_by_id[1]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [0, 17, 19, 25, 1, 18, 20, 26, 30, 69, 74, 32, 34, 71, 75, 36, 27, 68, 73, 80]
        """
        conn_field = self._meshed_region.elements.connectivities_field
        nodes_scoping = self._meshed_region.nodes.scoping
        return ConnectivityListByIndex(
            field=conn_field, mode=ReturnMode.IDX, scoping=nodes_scoping
        )

    @property
    def node_to_element_connectivity(self) -> ConnectivityListByIndex:
        """Returns a connectivity map between node index and element indexes.

        To get the connectivity map by node ID, use the 'by_id' property of the object returned.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> node_to_element_connectivity = simulation.mesh.node_to_element_connectivity
        >>> print(node_to_element_connectivity[0]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [0, 1, 2, 3, 4, 5, 6, 7]
        >>> node_to_element_connectivity_by_id = node_to_element_connectivity.by_id
        >>> print(node_to_element_connectivity_by_id[1]) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [0, 1, 2, 3, 4, 5, 6, 7]
        """
        conn_field = self._meshed_region.nodes.nodal_connectivity_field
        elems_scoping = self._meshed_region.elements.scoping
        return ConnectivityListByIndex(
            field=conn_field, mode=ReturnMode.IDX, scoping=elems_scoping
        )

    @property
    def face_to_node_ids_connectivity(self) -> ConnectivityListByIndex:
        """Returns a connectivity map between face index and node IDs.

        To get the connectivity map by face ID, use the 'by_id' property of the object returned.
        """
        conn_field = self._meshed_region.faces.faces_nodes_connectivity_field
        nodes_scoping = self._meshed_region.nodes.scoping
        return ConnectivityListByIndex(
            field=conn_field, mode=ReturnMode.IDS, scoping=nodes_scoping
        )

    @property
    def face_to_node_connectivity(self) -> ConnectivityListByIndex:
        """Returns a connectivity map between face index and node indexes.

        To get the connectivity map by face ID, use the 'by_id' property of the object returned.
        """
        conn_field = self._meshed_region.faces.faces_nodes_connectivity_field
        nodes_scoping = self._meshed_region.nodes.scoping
        return ConnectivityListByIndex(
            field=conn_field, mode=ReturnMode.IDX, scoping=nodes_scoping
        )

    @property
    def unit(self) -> str:
        """Returns the unit of the mesh (same as coordinates of the mesh).

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.unit) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        m
        """
        return self._meshed_region.unit

    @unit.setter
    def unit(self, value: str):
        """Set the unit of the mesh.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh.unit) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        m
        >>> simulation.mesh.unit = "mm"
        >>> print(simulation.mesh.unit) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        mm
        """
        self._meshed_region.unit = value

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.MeshedRegion` object."""
        return self._meshed_region

    def plot(self, **kwargs):
        """Plots the Mesh.

        Parameters
        ----------
        kwargs:
            Additional keyword arguments for the plotter. For additional keyword
            arguments, see ``help(pyvista.plot)``.

        Returns
        -------
        A Plotter instance of the current plotting back-end.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> from ansys.dpf.post.common import elemental_properties
        >>> example_path = examples.download_all_kinds_of_complexity()
        >>> simulation = post.StaticMechanicalSimulation(example_path)
        >>> mesh = simulation.mesh
        >>> mesh.plot()

        """
        return self._core_object.plot(**kwargs)

    @property
    def coordinates(self) -> post.DataFrame:
        """Returns the nodal coordinates of the Mesh as a DataFrame.

        Examples
        --------
        >>> from ansys.dpf.post import StaticMechanicalSimulation
        >>> from ansys.dpf.post import examples
        >>> simulation = StaticMechanicalSimulation(examples.static_rst)
        >>> mesh = simulation.mesh
        >>> coord = mesh.coordinates
        >>> print(coord)  # doctest: +NORMALIZE_WHITESPACE
                     results  coord (m)
         node_ids components
                1          X 1.5000e-02
                           Y 4.5000e-02
                           Z 1.5000e-02
                2          X 1.5000e-02
                           Y 4.5000e-02
                           Z 0.0000e+00
              ...        ...        ...
        """
        from ansys.dpf.post.simulation import vector_component_names

        label = "coord"
        fields_container = dpf.FieldsContainer()
        fields_container.add_field(
            label_space={}, field=self._core_object.nodes.coordinates_field
        )
        return post.DataFrame(
            data=fields_container,
            index=index.MultiIndex(
                indexes=[
                    index.MeshIndex(
                        location=locations.nodal,
                        scoping=self._core_object.nodes.scoping,
                        fc=fields_container,
                    ),
                    index.CompIndex(values=vector_component_names),
                ]
            ),
            columns=index.MultiIndex(
                indexes=[
                    index.ResultsIndex(values=[label], units=fields_container[0].unit)
                ]
            ),
        )
