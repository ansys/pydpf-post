"""This module contains ElementList, ElementType and Element classes."""

from __future__ import annotations

from collections.abc import Collection, Iterator
from typing import List, Union

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

    def __str__(self):
        """Returns a string representation of the Element Type."""
        return self.name

    def __repr__(self):
        """Returns a string representation of the Element Type."""
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
    def num_nodes(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.n_nodes`."""
        return self._resolve().n_nodes

    @property
    def type_info(self) -> ElementType:
        """Gets an element descriptor, See :py:meth:`ansys.dpf.core.elements.Element.id`."""
        return ElementType(self._resolve().type.value)

    @property
    def type(self) -> int:
        """Returns the ID of the Element Type."""
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
        """Returns string representation of an Element."""
        return self._resolve().__repr__()

    def __str__(self) -> str:
        """Returns string representation of an Element."""
        return self._resolve().__str__()


class ElementListIterator(Iterator):
    """Iterator class for the ElementList."""

    def __init__(self, el_list: ElementListIdx):
        """Constructs an Iterator from an element list."""
        self._el_list = el_list
        self._idx = 0

    def __next__(self) -> Element:
        """Returns the next Element in the list."""
        if self._idx >= self._el_list.__len__():
            raise StopIteration

        ret = self._el_list[self._idx]
        self._idx += 1
        return ret

    def __iter__(self) -> Iterator:
        """Returns a new Iterator object."""
        return ElementListIterator(self._el_list)


class ElementListIdx(Collection):
    """List of Elements."""

    def __init__(self, elements: elements.Elements):
        """Constructs list from existing dpf.core.elements.Elements list."""
        self._elements = elements

    def __getitem__(self, idx: int) -> Element:
        """Delegates to element_by_id() if by_id, otherwise to element_by_index()."""
        return Element(self._elements, idx)

    def __contains__(self, el: Element) -> bool:
        """Checks if the given element in the list."""
        return el.index >= 0 and el.index < self.__len__()

    def __iter__(self) -> ElementListIterator:
        """Returns an Iterator object on the list."""
        return ElementListIterator(self)

    def __len__(self) -> int:
        """Returns the number of elements in the list."""
        return self._elements.n_elements

    @property
    def by_id(self) -> ElementListById:
        """Returns an equivalent list accessible with ID instead of index."""
        return ElementListById(self._elements)


class ElementListById(ElementListIdx):
    """Wrapper class for accessing Elements by ID instead of index."""

    def __init__(self, elements: elements.Elements):
        """Constructs an ElementListById from an Elements instance."""
        super().__init__(elements)

    def __getitem__(self, id: int) -> Element:
        """Access an Element with an ID."""
        idx = self._elements.scoping.index(id)
        return super().__getitem__(idx)

    def __contains__(self, el: Element):
        """Checks if the given element is in the list."""
        return el.id in self._elements.scoping.ids

    def __iter__(self) -> ElementListIterator:
        """Returns an iterator object on the list."""
        return super().__iter__()

    def __len__(self) -> int:
        """Returns the number of elements in the list."""
        return self._elements.n_elements
