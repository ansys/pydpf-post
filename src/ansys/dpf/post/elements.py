"""This module contains ElementList, ElementType and Element classes."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Dict, List, Union

import ansys.dpf.core as dpf
import ansys.dpf.core.elements as elements
import ansys.dpf.core.nodes as nodes  # noqa: F401

import ansys.dpf.post as post
from ansys.dpf.post import index, locations
from  ansys.dpf.post.fields_container import PropertyFieldsContainer

class ElementType(dpf.ElementDescriptor):
    """Wrapper type to instantiate an ElementDescriptor from an int."""

    def __init__(self, arg: Union[dpf.ElementDescriptor, int]):
        """Constructs an ElementType from an existing descriptor or enum id."""
        _obj = arg
        if isinstance(arg, int):
            _obj = dpf.element_types.descriptor(arg)

        if not isinstance(_obj, dpf.ElementDescriptor):
            raise TypeError(f"Given argument is not an int nor an ElementDescriptor")

        super().__init__(
            _obj.enum_id,
            _obj.description,
            _obj.name,
            _obj.shape,
            _obj.n_corner_nodes,
            _obj.n_mid_nodes,
            _obj.n_nodes,
            _obj.is_solid,
            _obj.is_shell,
            _obj.is_beam,
            _obj.is_quadratic,
        )

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.__str__()

class Element:
    """Proxy class wrapping dpf.core.elements.Element."""

    def __init__(self, elements: elements.Elements, index: int):
        """Constructs a Proxy Element object."""
        self._elements = elements
        self._index = index

    def _resolve(self):
        """Returns the original Element object in the original list."""
        return self._elements[self._index]

    @property
    def node_ids(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.elements.Element.node_ids`."""
        return self._resolve().node_ids

    @property
    def id(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.id`."""
        return self._resolve().id

    @property
    def index(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.index`."""
        return self._resolve().index

    @property
    def nodes(self) -> List[nodes.Node]:
        """See :py:meth:`ansys.dpf.core.elements.Element.nodes`."""
        return self._resolve().nodes

    @property
    def n_nodes(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.n_nodes`."""
        return self._resolve().n_nodes

    @property
    def type_info(self) -> ElementType:
        """Gets an element descriptor, See :py:meth:`ansys.dpf.core.elements.Element.id`."""
        return ElementType(self._resolve().type.value)

    @property
    def type_id(self) -> int:
        return self._resolve().type.value

    @property
    def shape(self) -> str:
        """See :py:meth:`ansys.dpf.core.elements.Element.shape`."""
        return self._resolve().shape

    @property
    def connectivity(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.elements.Element.connectivity`."""
        return self._resolve().connectivity
    
    def __repr__(self) -> str:
        return self._resolve().__repr__()
    
    def __str__(self) -> str:
        return self._resolve().__str__()

class ElementListIdx(Sequence):
    """List of Elements."""

    def __init__(self, elements: elements.Elements):
        """Constructs list from existing dpf.core.elements.Elements list."""
        self._elements = elements

    def __getitem__(self, idx: int) -> Element:
        """Delegates to element_by_id() if by_id, otherwise to element_by_index()."""
        return Element(self._elements, idx)

    def __len__(self) -> int:
        """Returns the number of elements in the list."""
        return self._elements.n_elements

    @property
    def by_id(self) -> ElementListById:
        return ElementListById(self._elements)

    @property
    def types(self) -> post.DataFrame:
        """Returns mapping of element id to corresponding type."""
        field: dpf.Field = self._elements.element_types_field
        label = "el_type_id"
        fields_container = PropertyFieldsContainer()
        fields_container.add_field(
            label_space={}, field=field
        )

        return post.DataFrame(
            data=fields_container,
            index=index.MultiIndex(
                indexes=[
                    index.MeshIndex(
                        location=locations.elemental,
                        scoping=self._elements.scoping,
                        fc=fields_container
                    )
                ]
            ),
            columns=index.MultiIndex(
                indexes=[
                    index.ResultsIndex(values=[label])
                ]
            )
        )

class ElementListById(ElementListIdx):
    def __init__(self, elements: elements.Elements):
        super().__init__(elements)

    def __getitem__(self, id: int) -> Element:
        idx = self._elements.scoping.index(id)
        return super().__getitem__(idx)
