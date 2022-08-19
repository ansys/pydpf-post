"""This module contains a class describing the result objects."""

from ansys.dpf.core.common import locations

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.post.dpf_path import DpfPath


class Definition:
    """Class containing the attributes as properties of the result object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        self._location = locations.nodal
        self._element_scoping = None
        self._node_scoping = None
        self._named_selection = None
        self._grouping = None
        self._mapdl_grouping = None
        self._time_scoping = None
        self._time = None
        self._set = None
        self._path = None

        self.__location_locked = False
        self.__element_scoping_locked = False

        if _AvailableKeywords.location in kwargs:
            self._location = kwargs.pop(_AvailableKeywords.location)
        if _AvailableKeywords.element_scoping in kwargs:
            self._element_scoping = kwargs.pop(_AvailableKeywords.element_scoping)
        if _AvailableKeywords.node_scoping in kwargs:
            self._node_scoping = kwargs.pop(_AvailableKeywords.node_scoping)
        if _AvailableKeywords.named_selection in kwargs:
            self._named_selection = kwargs.pop(_AvailableKeywords.named_selection)
        if _AvailableKeywords.grouping in kwargs:
            self._grouping = kwargs.pop(_AvailableKeywords.grouping)
        if _AvailableKeywords.mapdl_grouping in kwargs:
            self._mapdl_grouping = kwargs.pop(_AvailableKeywords.mapdl_grouping)
        if _AvailableKeywords.time_scoping in kwargs:
            self._time_scoping = kwargs.pop(_AvailableKeywords.time_scoping)
        if _AvailableKeywords.time in kwargs:
            self._time = kwargs.pop(_AvailableKeywords.time)
        if _AvailableKeywords.set in kwargs:
            self._set = kwargs.pop(_AvailableKeywords.set)
        if _AvailableKeywords.path in kwargs:
            self._path = kwargs.pop(_AvailableKeywords.path)

        if len(kwargs) > 0:
            txt = "Following keyword argument(s) could not be computed: \n"
            for key in kwargs:
                txt += "- " + key + "\n"
                txt += "\n"
                txt += "Use 'post.print_available_keywords()' to see which keywords can be used."
            raise Exception(txt)

    def __str__(self):
        """Return the string representation."""
        txt = "Object properties:\n"
        for attr, value in vars(self).items():
            if not (
                attr.startswith("__")
                or attr.startswith("_data_sources")
                or attr.startswith("_model")
                or attr.startswith("_Definition")
            ):
                if value is not None:
                    if attr.startswith("_"):
                        attr = attr[1:]
                    if isinstance(value, int):
                        value_temp = value
                        value = str(value_temp)
                    if isinstance(value, float):
                        value_temp = value
                        value = str(value_temp)
                    if isinstance(value, list):
                        value_temp = "["
                        for val in value:
                            value_temp += str(val) + ","
                        value_temp = value_temp[: len(value_temp) - 1]
                        value_temp += "]"
                        value = value_temp
                    txt += f" - {attr:10} : " + value + "\n"
        return txt

    @property
    def location(self) -> str:
        """Return or set the location.

        Accepts :class:`ansys.dpf.core.common.locations`.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.location = post.locations.elemental
        >>> stress.location
        'Elemental'

        """
        return self._location

    @location.setter
    def location(self, value):
        if self.__location_locked:
            raise Exception(
                "Location can not be set outside of the instantiation "
                "of the result object in this case."
            )
        if value is not None:
            if not isinstance(value, str):
                raise TypeError("Expected type is str.")
        self._location = value

    @property
    def element_scoping(self):
        """Return or set elemental scoping.

        Accepts sequence or :class:`ansys.dpf.core.Scoping`.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.element_scoping = [1, 4]
        """
        return self._element_scoping

    @element_scoping.setter
    def element_scoping(self, value):
        if self.__element_scoping_locked:
            raise Exception(
                "Element scoping can not be set outside of the "
                "instantiation of the result object in this case."
            )
        self._element_scoping = value

    @property
    def node_scoping(self):
        """Return or set nodal scoping property.

        Accepts sequencye or :class:`ansys.dpf.core.Scoping`.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.node_scoping = [1, 4]
        """
        return self._node_scoping

    @node_scoping.setter
    def node_scoping(self, value):
        self._node_scoping = value

    @property
    def named_selection(self) -> str:
        """Return or set the named selection.

        MAPDL named selections are in upper-case.  Lower case name selections
        will automatically be set to uppercase.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.node_scoping = 'SELECTION'
        """
        return self._named_selection

    @named_selection.setter
    def named_selection(self, value):
        if not isinstance(value, str):
            raise TypeError("Expected type is str.")
        self._named_selection = value.upper()

    @property
    def grouping(self) -> str:
        """Return or set the grouping.

        Accepts :class:`ansys.dpf.post.grouping` or str.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.node_scoping = post.grouping.by_el_shape
        """
        return self._grouping

    @grouping.setter
    def grouping(self, value):
        if not isinstance(value, str):
            raise TypeError("Expected type is str.")
        self._grouping = value

    @property
    def time_scoping(self):
        """Return or set the time scoping.

        Accepts sequence or :class:`ansys.dpf.core.Scoping`.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
        >>> stress = solution.stress()
        >>> stress.time_scoping = [1, 4]
        """
        return self._time_scoping

    @time_scoping.setter
    def time_scoping(self, value):
        self._time_scoping = value

    @property
    def mapdl_grouping(self) -> int:
        """Return or set the MAPDL grouping.

        Examples
        --------
        Set the grouping to the ``SOLID186`` MAPDL element type.

        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.mapdl_grouping = 186
        """
        return self._mapdl_grouping

    @mapdl_grouping.setter
    def mapdl_grouping(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected type is int.")
        self._mapdl_grouping = value

    @property
    def time(self) -> float:
        """Return or set the time step.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
        >>> stress = solution.stress()
        >>> stress.time = 0.2
        """
        return self._time

    @time.setter
    def time(self, value):
        if not isinstance(value, float):
            raise TypeError("Expected type is float.")
        self._time = value

    @property
    def set(self) -> int:
        """Return or set the result set.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> stress.set = 2
        """
        return self._set

    @set.setter
    def set(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected type is int.")
        self._set = value

    @property
    def path(self) -> DpfPath:
        """Return or set the coordinates path.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> stress = solution.stress()
        >>> coordinates = [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.0, 0.1, 0.0]]
        >>> path = post.create_path_on_coordinates(coordinates=coordinates)
        >>> stress.path = path
        """
        return self._path

    @path.setter
    def path(self, value):
        self._path = value
