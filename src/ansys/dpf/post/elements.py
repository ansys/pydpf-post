"""This module contains ElementList, ElementType and Element classes."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Dict, List, Union

import ansys.dpf.core as dpf
import ansys.dpf.core.elements as elements
import ansys.dpf.core.nodes as nodes  # noqa: F401


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
    def type(self) -> ElementType:
        """Gets an element descriptor, See :py:meth:`ansys.dpf.core.elements.Element.id`."""
        return ElementType(self._resolve().type.value)

    @property
    def shape(self) -> str:
        """See :py:meth:`ansys.dpf.core.elements.Element.shape`."""
        return self._resolve().shape

    @property
    def connectivity(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.elements.Element.connectivity`."""
        return self._resolve().connectivity


class ElementList(Sequence):
    """List of Elements."""

    def __init__(self, elements: elements.Elements, by_id=True):
        """Constructs list from existing dpf.core.elements.Elements list."""
        self._elements = elements
        self.by_id = by_id

    def __getitem__(self, key: int) -> Element:
        """Delegates to element_by_id() if by_id, otherwise to element_by_index()."""
        index = key
        if self.by_id:
            index = self._elements.element_by_id(key).index

        return Element(self._elements, index)

    def __len__(self) -> int:
        """Returns the number of elements in the list."""
        return self._elements.n_elements

    # def __repr__(self) -> str:
    #    return list(iter(self._meshed_region.elements)).__repr__()

    @property
    def types(self) -> Dict[int, ElementType]:
        """Returns mapping of element id to corresponding type."""
        # TODO: Dataframe
        field: dpf.Field = self._elements.element_types_field
        keys = field.scoping.ids

        int_to_ed = lambda i: ElementType(int(i))
        vals = map(int_to_ed, field.data)
        return dict(zip(keys, vals))
