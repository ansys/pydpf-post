"""Module containing the ``Mesh`` class.

Mesh
----

"""
from __future__ import annotations

from typing import List

import ansys.dpf.core as dpf
from ansys.dpf.core.nodes import Node

import ansys.dpf.post as post
from ansys.dpf.post import index, locations
from ansys.dpf.post.connectivity import ConnectivityListIdx, ReturnMode
from ansys.dpf.post.elements import Element, ElementListIdx
from ansys.dpf.post.named_selection import NamedSelectionsDict
from ansys.dpf.post.nodes import NodeListIdx
from ansys.dpf.post.prop_fields_container import _PropertyFieldsContainer


class Mesh:
    """Exposes the complete mesh of the simulation."""

    def __init__(self, meshed_region: dpf.MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    def __str__(self):
        """String representation of this class."""
        return str(self._meshed_region).replace("Meshed Region", "Mesh")

    @property
    def named_selections(self) -> NamedSelectionsDict:
        """Returns a mapping of scopings."""
        return NamedSelectionsDict(self._meshed_region)

    @property
    def node_ids(self) -> List[int]:
        """Returns the list of node IDs in the mesh."""
        return self._meshed_region.nodes.scoping.ids

    @property
    def num_nodes(self) -> int:
        """Returns the number of nodes in the mesh."""
        return len(self.node_ids)

    @property
    def element_ids(self) -> List[int]:
        """Returns the list of element IDs in the mesh."""
        return self._meshed_region.elements.scoping.ids

    @property
    def num_elements(self) -> int:
        """Returns the number of element in the mesh."""
        return len(self.element_ids)

    @property
    def elements(self) -> ElementListIdx:
        """Returns a list of elements indexed by ID."""
        return ElementListIdx(self._meshed_region.elements)

    def get_element_by_id(
        self, id: int  # pylint: disable=redefined-builtin
    ) -> Element:
        """Returns an element in the mesh from a given ID."""
        return self.elements.by_id[id]

    @property
    def nodes(self) -> NodeListIdx:
        """Returns a list of nodes indexed by ID."""
        return NodeListIdx(self._meshed_region.nodes)

    def get_node_by_id(self, id: int) -> Node:  # pylint: disable=redefined-builtin
        """Returns a node of the mesh from a given ID."""
        return self.nodes.by_id[id]

    @property
    def element_types(self) -> post.DataFrame:
        """Returns a DataFrame containing element types ID."""
        label = "elem_type_id"
        fields_container = _PropertyFieldsContainer()
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
        """Returns a DataFrame containing element materials ID."""
        label = "material_id"
        fields_container = _PropertyFieldsContainer()
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
    def element_to_node_ids_connectivity(self):
        """Returns a."""
        conn_field = self._meshed_region.elements.connectivities_field
        nodes_scoping = self._meshed_region.nodes.scoping
        return ConnectivityListIdx(
            field=conn_field, mode=ReturnMode.IDS, scoping=nodes_scoping
        )

    @property
    def node_to_element_ids_connectivity(self):
        """Returns a list of Element ID for a given Node index."""
        conn_field = self._meshed_region.nodes.nodal_connectivity_field
        elems_scoping = self._meshed_region.elements.scoping
        return ConnectivityListIdx(
            field=conn_field, mode=ReturnMode.IDS, scoping=elems_scoping
        )

    @property
    def element_to_node_connectivity(self):
        """Returns a list of Node index for a given Element index."""
        conn_field = self._meshed_region.elements.connectivities_field
        nodes_scoping = self._meshed_region.nodes.scoping
        return ConnectivityListIdx(
            field=conn_field, mode=ReturnMode.IDX, scoping=nodes_scoping
        )

    @property
    def node_to_element_connectivity(self):
        """Returns a list of Element index for a given Node index."""
        conn_field = self._meshed_region.nodes.nodal_connectivity_field
        elems_scoping = self._meshed_region.elements.scoping
        return ConnectivityListIdx(
            field=conn_field, mode=ReturnMode.IDX, scoping=elems_scoping
        )

    @property
    def unit(self) -> str:
        """Returns the unit of the mesh (same as coordinates of the mesh)."""
        return self._meshed_region.unit

    @unit.setter
    def unit(self, value: str):
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
                     results   coord (m)
        node_ids  components
               1           X  1.5000e-02
                           Y  4.5000e-02
                           Z  1.5000e-02
               2           X  1.5000e-02
                           Y  4.5000e-02
                           Z  0.0000e+00
             ...

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
