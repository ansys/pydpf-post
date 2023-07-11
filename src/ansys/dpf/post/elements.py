"""This module contains ElementList, ElementType and Element classes."""

from __future__ import annotations

from collections.abc import Collection, Iterator
from typing import List, Union

import ansys.dpf.core as dpf
from ansys.dpf.core import elements
from ansys.dpf.core.nodes import Node


class ElementType:
    """Wrapper type to instantiate an ElementDescriptor from an int."""

    def __init__(self, arg: Union[dpf.ElementDescriptor, int]):
        """Constructs an ElementType from an existing descriptor or enum id."""
        self._el_desc = arg
        if isinstance(arg, int):
            self._el_desc = dpf.element_types.descriptor(arg)

        if not isinstance(self._el_desc, dpf.ElementDescriptor):
            raise TypeError("Given argument is not an int nor an ElementDescriptor")

    @property
    def enum_id(self) -> dpf.element_types:
        """Element type in the element_types enum."""
        return self._el_desc.enum_id

    @property
    def description(self) -> str:
        """Specifies the element geometry and integration order."""
        return self._el_desc.description

    @property
    def name(self) -> str:
        """Short name of the element type."""
        return self._el_desc.name

    @property
    def shape(self) -> str:
        """Can be ``"solid"``,``"shell"`` or ``"beam"``."""
        return self._el_desc.shape

    @property
    def num_corner_nodes(self) -> int:
        """Returns the number of corner nodes."""
        return self._el_desc.n_corner_nodes

    @property
    def num_mid_nodes(self) -> int:
        """Returns the number of middle nodes."""
        return self._el_desc.n_mid_nodes

    @property
    def num_nodes(self) -> int:
        """Returns the total number of nodes."""
        return self._el_desc.n_nodes

    @property
    def is_solid(self) -> bool:
        """Whether the element is a solid."""
        return self._el_desc.is_solid

    @property
    def is_shell(self) -> bool:
        """Whether the element is a shell."""
        return self._el_desc.is_shell

    @property
    def is_beam(self) -> bool:
        """Whether the element is a beam."""
        return self._el_desc.is_beam

    @property
    def is_quadratic(self) -> bool:
        """Whether the element is quadratic."""
        return self._el_desc.is_quadratic

    def __str__(self) -> str:
        """Returns a string representation of the Element Type."""
        return self._el_desc.__str__().replace(
            "Element descriptor\n------------------", "Element Type\n------------"
        )

    def __repr__(self):
        """Returns a string representation of the Element Type."""
        return self.__str__()


class Element:
    """Proxy class wrapping dpf.core.elements.Element."""

    def __init__(self, element: elements.Element):
        """Constructs a Proxy Element object."""
        self._element = element

    @property
    def node_ids(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.elements.Element.node_ids`."""
        return self._element.node_ids

    @property
    def id(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.id`."""
        return self._element.id

    @property
    def index(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.index`."""
        return self._element.index

    @property
    def nodes(self) -> List[Node]:
        """See :py:meth:`ansys.dpf.core.elements.Element.nodes`."""
        return self._element.nodes

    @property
    def num_nodes(self) -> int:
        """See :py:meth:`ansys.dpf.core.elements.Element.n_nodes`."""
        return self._element.n_nodes

    @property
    def type_info(self) -> ElementType:
        """Gets an element descriptor, See :py:meth:`ansys.dpf.core.elements.Element.id`."""
        return ElementType(self._element.type.value)

    @property
    def type(self) -> elements.element_types:
        """Returns the Element Type."""
        return self._element.type

    @property
    def shape(self) -> str:
        """See :py:meth:`ansys.dpf.core.elements.Element.shape`."""
        return self._element.shape

    @property
    def to_node_connectivity(self) -> List[int]:
        """See :py:meth:`ansys.dpf.core.elements.Element.connectivity`."""
        return self._element.connectivity

    def __repr__(self) -> str:
        """Returns string representation of an Element."""
        return f"Element(type={self.type},index={self.index},id={self.id},shape={self.shape})"

    def __str__(self) -> str:
        """Returns string representation of an Element."""
        return self._element.__str__()


class ElementListIterator(Iterator):
    """Iterator class for the ElementList."""

    def __init__(self, el_list: elements.Elements):
        """Constructs an Iterator from an element list."""
        self._el_list = el_list
        self._idx = 0

    def __next__(self) -> Element:
        """Returns the next Element in the list."""
        if self._idx >= self._el_list.__len__():
            raise StopIteration

        ret = Element(self._el_list, self._idx)
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

    def __getitem__(self, index: int) -> Element:
        """Delegates to element_by_id() if by_id, otherwise to element_by_index()."""
        return Element(self._elements, index)

    def __contains__(self, el: Element) -> bool:
        """Checks if the given element in the list."""
        return el.index >= 0 and el.index < self.__len__()

    def __iter__(self) -> ElementListIterator:
        """Returns an Iterator object on the list."""
        return ElementListIterator(self._elements)

    def __len__(self) -> int:
        """Returns the number of elements in the list."""
        return self._elements.n_elements

    @property
    def by_id(self) -> ElementListById:
        """Returns an equivalent list accessible with ID instead of index."""
        return ElementListById(self._elements)

    def _short_list(self) -> str:
        _str = "["
        if self.__len__() > 3:
            _fst = Element(self._elements, 0).type_info.name
            _lst = Element(self._elements, self.__len__() - 1).type_info.name
            _str += f"{_fst}, ..., {_lst}"
        else:
            el_list = [Element(self._elements, idx) for idx in range(self.__len__())]
            _str += ", ".join(map(lambda el: el.type_info.name, el_list))
        _str += "]"
        return _str

    def __str__(self) -> str:
        """Returns a string representation of ElementListIdx."""
        return self._short_list()

    def __repr__(self) -> str:
        """Returns a string representation of ElementListIdx."""
        return f"ElementListIdx({self.__str__()}, __len__={self.__len__()})"


class ElementListById(ElementListIdx):
    """Wrapper class for accessing Elements by ID instead of index."""

    def __init__(self, elements: elements.Elements):
        """Constructs an ElementListById from an Elements instance."""
        super().__init__(elements)

    def __getitem__(self, id: int) -> Element:  # pylint: disable=redefined-builtin
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

    def __repr__(self):
        """Returns a string representation of ElementListById."""
        return f"ElementListById({super().__str__()}, __len__={self.__len__()})"
