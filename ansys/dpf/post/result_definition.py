"""This module contains a class describing the result objects."""

from ansys.dpf.core.common import locations
from ansys.dpf.post.common import _AvailableKeywords

class Definition:
    """Class containing the attributes as property (with setter and getter)
    of the result object."""

    def __init__(self, **kwargs):
        self._location = locations.nodal
        self._element_scoping = None
        self._node_scoping = None
        self._named_selection = None
        self._grouping = None
        self._mapdl_grouping = None
        self._time_scoping = None
        self._time = None
        self._set = None

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

        if (len(kwargs) > 0):
            txt = "Following keyword argument(s) could not be computed: \n"
            for key in kwargs:
                txt += "- " + key + "\n"
                txt += "\n"
                txt += "Use 'post.print_available_keywords()' to see which keywords can be used."
            raise Exception(txt)

    def __str__(self):
        txt = "Object properties:\n"
        for attr, value in vars(self).items():
            if not (attr.startswith('__') or attr.startswith('_data_sources') or attr.startswith('_model') or attr.startswith('_Definition')):
                if value is not None:
                    if attr.startswith('_'):
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
                        value_temp = value_temp[:len(value_temp)-1]
                        value_temp += "]"
                        value = value_temp
                    txt += f" - {attr:10} : " + value + "\n"
        return txt

    @property
    def location(self) -> str:
        """str: Location property. Can be set.
        post.locations enum can be used.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
        >>> stress = solution.stress()
        >>> stress.location = post.locations.elemental
        """
        return self._location

    @location.setter
    def location(self, value):
        if self.__location_locked:
            raise Exception("Location can not be set outside of the instantiation of the result object in this case.")
        if value is not None:
            if not isinstance(value, str):
                raise TypeError("Expected type is str.")
        self._location = value

    @property
    def element_scoping(self):
        """Elemental scoping property. Can be set.
        Available types: list, or dpf.core.Scoping.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
        >>> stress = solution.stress()
        >>> stress.element_scoping = [1, 4]
        """
        return self._element_scoping

    @element_scoping.setter
    def element_scoping(self, value):
        if self.__element_scoping_locked:
            raise Exception("Element scoping can not be set outside of the instantiation of the result object in this case.")
        self._element_scoping = value

    @property
    def node_scoping(self):
        """Nodal scoping property. Can be set.
        Available types: list, or dpf.core.Scoping.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
        >>> stress = solution.stress()
        >>> stress.node_scoping = [1, 4]
        """
        return self._node_scoping

    @node_scoping.setter
    def node_scoping(self, value):
        self._node_scoping = value

    @property
    def named_selection(self) -> str:
        """str: Named selection property. Can be set.

        MAPDL named selections are in upper-case.  Lower case name
        selections will automatically be set to uppercase.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
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
        """str: Grouping property. Can be set.
        post.grouping enum can be used.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
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
        """Time scoping property. Can be set.
        Available types: list, or dpf.core.Scoping.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
        >>> stress = solution.stress()
        >>> stress.time_scoping = [1, 4]
        """
        return self._time_scoping

    @time_scoping.setter
    def time_scoping(self, value):
        self._time_scoping = value

    @property
    def mapdl_grouping(self) -> int:
        """int: Mapdl grouping property. Can be set.
        Following is an example to get only solid 186
        mapdl element type.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
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
        """float: time step property. Can be set.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
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
        """int: Set property. Can be set.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> solution = post.solution("file.rst")
        >>> stress = solution.stress()
        >>> stress.set = 2
        """
        return self._set

    @set.setter
    def set(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected type is int.")
        self._set = value
