# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
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

"""Module containing the ``DpfPath`` class."""
from ansys.dpf.core import Field
from ansys.dpf.core.common import locations, natures
import numpy as np


def create_path_on_coordinates(coordinates):
    """Create a DPF path object.

    You can use this path object to request results on a specific path of
    coordinates.

    Parameters
    ----------
    coordinates : list[list[int]], field, numpy.ndarray
        3D coordinates.

    Examples
    --------
    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> coordinates = [[0.024, 0.03, 0.003]]
    >>> for i in range(1, 51):
    ...     coord_copy = coordinates[-1].copy()
    ...     coord_copy[1] = coord_copy[0] + i * 0.001
    ...     coordinates.append(coord_copy)
    >>> path_on_coord = post.create_path_on_coordinates(
    ... coordinates=coordinates
    ... )
    >>> solution = post.load_solution(examples.static_rst)
    >>> stress = solution.stress(path=path_on_coord)

    """
    return DpfPath(coordinates=coordinates)


class DpfPath:
    """Describes a set of coordinates.

    Parameters
    ----------
    coordinates : list[list[int]], field, arrays
        3D coordinates.

    Example
    -------
    Create coordinates from a list.

    >>> from ansys.dpf import post
    >>> coordinates = [[0.024, 0.03, 0.003]]
    >>> for i in range(1, 51):
    ...     coord_copy = coordinates[-1].copy()
    ...     coord_copy[1] = coord_copy[0] + i * 0.001
    ...     coordinates.append(coord_copy)
    >>> dpf_path = post.create_path_on_coordinates(coordinates=coordinates)

    Create coordinates from a :class:`numpy.ndarray`.

    >>> import numpy as np
    >>> coordinates = np.empty((50, 3) )
    >>> coordinates[:] = [0.024, 0.03, 0.003]
    >>> coordinates[:, 1] = np.linspace(0.03, 0.74)
    >>> dpf_path = post.create_path_on_coordinates(coordinates=coordinates)

    """

    def __init__(self, coordinates):
        """Initialize this class."""
        if isinstance(coordinates, Field):
            self._field = coordinates
        else:
            coord_length = len(coordinates)
            if isinstance(coordinates, list):
                if isinstance(coordinates[0], float):
                    coord_length /= 3
            elif isinstance(coordinates, (np.ndarray, np.generic)):
                if len(coordinates.shape) == 1:
                    coord_length /= 3
            self._field = Field(nature=natures.vector, location=locations.nodal)
            self._field.scoping.ids = list(range(1, int(coord_length) + 1))
            self._field.data = coordinates

    @property
    def coordinates(self):
        """Coordinates of the path."""
        return self._field.data
