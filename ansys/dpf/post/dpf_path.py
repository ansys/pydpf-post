"""Module containing the DpfPath class that will
define a path of coordinates to set the result on.
"""

from ansys.dpf.core import locations
from ansys.dpf.core import Scoping, Field

def create_path_on_coordinates(coordinates, scoping=None):
    """
    Create a dpf path object that can be used to request
    results on a specific path of coordinates.

    Parameters
    ----------
    coordinates : list[list[int]]
        List of coordinates.
    scoping : list[int]
        List of ids. Default is ``None``.
        In the case scoping is ``None``,
        scoping is set to
        list(range(1, len(coordinates) + 1)).
        Scoping location is set to nodal.

    Example
    -------
    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> coordinates = [[0.024, 0.03, 0.003]]
    >>> for i in range(1, 51):
    ...     coord_copy = ref.copy()
    ...     coord_copy[1] = coord_copy[0] + i * 0.001
    ...     coordinates.append(coord_copy)
    >>> path_on_coord = post.create_path_on_coordinates(
    ... coordinates=coordinates
    ... )
    >>> solution = post.load_solution(examples.static_rst)
    >>> stress = solution.stress(path=dpf_path)

    """
    return DpfPath(coordinates=coordinates, scoping=scoping)

class DpfPath:
    """This object describe a set of coordinates.
    It can be associated to a scoping (list of ids)."""

    def __init__(self, coordinates, scoping=None):
        """
        DpfPath object constructor.

        Parameters
        ----------
        coordinates : list[list[int]]
            List of coordinates.
        scoping : list[int]
            List of ids. Default is ``None``.
            In the case scoping is ``None``,
            scoping is set to
            list(range(1, len(coordinates) + 1)).
            Scoping location is set to nodal.

        Example
        -------
        >>> coordinates = [[0.024, 0.03, 0.003]]
        >>> for i in range(1, 51):
        ...     coord_copy = ref.copy()
        ...     coord_copy[1] = coord_copy[0] + i * 0.001
        ...     coordinates.append(coord_copy)
        >>> dpf_path = post.DpfPath(coordinates=coordinates)

        """
        # check inputs
        error_text = """coordinates attribute must
        be a list of list of int."""
        if not isinstance(coordinates, list):
            raise ValueError(error_text)
        else:
            if len(coordinates) > 0 and not isinstance(coordinates[0], list):
                raise ValueError(error_text)
        _scoping = None
        if scoping is None:
            _scoping = Scoping(location=locations.nodal)
            _scoping.ids = list(range(1, len(coordinates) + 1))
        else:
            error_text = """scoping must be a list of int.
            Its size must be identical to the coordinates array one."""
            if isinstance(scoping, list) and len(scoping) == len(coordinates):
                if isinstance(scoping[0], int):
                    _scoping = Scoping(location=locations.nodal)
                    _scoping.ids = scoping
                else:
                    raise ValueError(error_text)
            elif isinstance(scoping, Scoping):
                _scoping = scoping
            else:
                raise ValueError(error_text)
        # create field
        self._field = Field(location=_scoping.location)
        self._field.data = coordinates
        self._field.scoping = _scoping

    @property
    def coordinates(self):
        return self._field.data

    @property
    def scoping(self):
        return self._field.scoping.ids