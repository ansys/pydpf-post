"""Module containing the _ResultData class ("result object") that
user will be able to use to compute through the DPF Post API.

This is a fields container wrapper."""

from textwrap import wrap

from ansys.dpf.core import FieldsContainer, Operator
from ansys.dpf.core.common import types, DefinitionLabels
from ansys.dpf.core.plotter import Plotter as DpfPlotter
from ansys.dpf.post.result_evaluation import ResultEvaluator
from ansys.dpf.post import errors as dpf_errors
from ansys.dpf.core import errors as core_errors


class ResultData:
    """Instance of the result of a dpf result object.
    Created thanks to the dpf result object.

    Parameters
    ----------
    Following list of keywords:
        - location
        - node_scoping
        - element_scoping
        - time
        - grouping

    The entire list of parameters can be found using
     ``post.print_available_keywords()``.

    Examples
    --------
    >>> from ansys.dpf import post
    >>> from ansys.dpf.post import examples
    >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
    >>> disp = solution.displacement()
    >>> disp_on_nodes = solution.displacement(node_scoping = [1, 23])
    >>> disp_on_named_selection = solution.displacement(named_selection="SELECTION")
    """

    def __init__(
        self,
        operator_name: str,
        data_sources,
        model,
        elem_average: bool,
        location: str = None,
        element_scoping=None,
        node_scoping=None,
        named_selection=None,
        time=None,
        grouping=None,
        phase=None,
        subresult=None,
        mapdl_grouping=None,
        set=None,
        path=None,
        time_scoping=None,
    ):

        self._evaluator = ResultEvaluator(
            operator_name,
            data_sources,
            model,
            elem_average,
            location,
            element_scoping,
            node_scoping,
            named_selection,
            time,
            grouping,
            phase,
            subresult,
            mapdl_grouping,
            set,
            path,
            time_scoping,
        )
        self._result_fields_container = None

    def __str__(self):
        self._evaluate_result()
        name = self.result_fields_container[0].name.split("_")
        txt = "%s result.\n" % name[0].capitalize()
        if self._evaluator.subresult is not None:
            txt += "%s component. \n" % self._evaluator.subresult
        txt += "\n"
        desc = "This result has been computed using dpf.core.Operator objects, which "
        desc += "have been chained together according to the following list:"
        txt += "\n".join(wrap(desc)) + "\n"
        for key, val in self._evaluator._chained_operators.items():
            txt += "- %s: " % key
            txt += val
            txt += "\n"
        # txt += "\n\n"
        # txt += self._evaluator._model.__str__()
        return txt

    def _evaluate_result(self):
        """First evaluation of the result."""
        if self._result_fields_container is None:
            self._result_fields_container = self._evaluator.evaluate_result()

    def _evaluate_result_forced(self):
        """Re-evaluation of the result."""
        self.result_fields_container = self._evaluator.evaluate_result()

    def get_all_label_spaces(self):
        """Returns all the label spaces contained in a result
        as a string.
        Labels can be used to select fields to plot.

        Returns
        -------
        list
            List of dictionary (list of label space)
        """
        self._evaluate_result()
        i = 0
        list_labels = []
        while i < len(self.result_fields_container):
            list_labels.append(self.result_fields_container.get_label_space(i))
            i += 1
        return list_labels

    @property
    def num_fields(self):
        """Returns the number of fields contained in the result."""
        self._evaluate_result()
        return len(self.result_fields_container)

    def get_data_at_field(self, field_index: int = 0):
        """Returns the data at the field with the specified index."""
        self._evaluate_result()
        # Needed to hold onto the field as a quick fix for a memory leak
        # which causes InProcess mode of PyDPF-Core 0.5.2 to crash
        # Requires a redesign
        owning_field = self.result_fields_container[field_index]
        data = owning_field.data
        try:
            data._owning_field = owning_field
        except AttributeError:
            pass
        return data

    def __getitem__(self, field_index: int = 0):
        """Override of the result item getter. Implements the fields_container_result
        item getter. Return the field at the field_index position in the result fields
        container.
        """
        self._evaluate_result()
        return self.result_fields_container[field_index]

    def get_scoping_at_field(self, field_index: int = 0):
        """Returns the scoping of the result."""
        self._evaluate_result()
        field = self.result_fields_container[field_index]
        return field.scoping.ids

    def _min_max(self, pin):
        """ "Chains the operators to compute the min/max values."""
        self._evaluate_result()
        max_operator = Operator("min_max_fc")
        max_operator.inputs.connect(self.result_fields_container)
        result = max_operator.get_output(pin, types.field)
        return result

    @property
    def max(self):
        """Returns the maximum values field.
        Chains the result operator to the "min_max_fc" operator,
        returns its result (output from pin 1).
        """
        return self._min_max(1)

    @property
    def max_data(self):
        """Returns the maximum values field data.
        Chains the result operator to the "min_max_fc" operator,
        returns its result (output from pin 1).
        """
        return self._min_max(1).data

    def get_max_data_at_field(self, field_index: int = 0):
        """Returns the maximum values field data at field_index.
        Chains the result operator to the "min_max_fc" operator,
        returns its result (output from pin 1).
        """
        return self._min_max(1).data[field_index]

    @property
    def min(self):
        """Returns the minimum values field.
        Chains the result operator to the "min_max_fc" operator,
        returns its result (output from pin 0).
        """
        return self._min_max(0)

    @property
    def min_data(self):
        """Returns the minimum values field data.
        Chains the result operator to the "min_max_fc" operator,
        returns its result (output from pin 0).
        """
        return self._min_max(0).data

    def get_min_data_at_field(self, field_index: int = 0):
        """Returns the minimum values field data at field_index.
        Chains the result operator to the "min_max_fc" operator,
        returns its result (output from pin 0).
        """
        return self._min_max(0).data[field_index]

    @property
    def result_fields_container(self):
        self._evaluate_result()
        return self._result_fields_container

    def _plot_contour_with_vtk_file(self):
        """Plot the contour result on its mesh support. The obtained figure depends on the
        support (can be a meshed_region or a time_freq_support).
        If transient analysis, plot the last result.

        This method is private, publishes a vtk file and print (using pyvista) from this file."""
        self._evaluate_result()
        pl = DpfPlotter(self._evaluator._model.metadata.meshed_region)
        pl._plot_contour_using_vtk_file(self.result_fields_container)

    def _sort_fields_container_with_labels(self, option_id, display_option):
        ids = option_id
        if isinstance(option_id, int):
            ids = [option_id]
        new_fields_container = FieldsContainer()
        label_spaces = []
        for opt_id in ids:
            i = 0
            while i < len(self.result_fields_container):
                label_space = {}
                label_space[display_option] = opt_id
                fc_label = self.result_fields_container.get_label_space(i)
                for lab in self.result_fields_container.labels:
                    if lab != display_option:
                        if display_option in fc_label:
                            if fc_label[display_option] == opt_id:
                                label_space[lab] = fc_label[lab]
                if len(fc_label) == len(label_space):
                    label_spaces.append(label_space)
                i += 1
        lab_names = []
        for lab_sp in label_spaces:
            for key, val in lab_sp.items():
                if key not in lab_names:
                    lab_names.append(key)
        for name in lab_names:
            new_fields_container.add_label(name)
        for label in label_spaces:
            field = self.result_fields_container.get_field(label)
            if not field:
                txt = """Arguments display_option/option_id are not correct,
                no corresponding field found to plot."""
                raise Exception(txt)
            new_fields_container.add_field(label, field)
        return new_fields_container

    def plot_contour(
        self,
        display_option: str = "time",
        option_id=1,
        **kwargs
    ):
        """Plot the contour result on its mesh support.

        The obtained figure depends on the support (can be a
        meshed_region or a time_freq_support).  If transient analysis,
        plot the last result if no time_scoping has been specified.
        The self.get_all_label_spaces() will return a string
        containing all the label spaces.

        Parameters
        ----------
        display_option : str, optional
            The name of the label you want to display. Default is ``"time"``.
        option_id: int, optional
             The list of label ids you want to display. Default is ``[1]``.
        **kwargs : optional
            Additional keyword arguments for the plotter.  See
            ``help(pyvista.plot)`` for additional keyword arguments.

        Examples
        --------
        Plot a result at the time_step number 1

        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.download_all_kinds_of_complexity())
        >>> stress = solution.stress(location=post.locations.nodal)
        >>> sx = stress.xx
        >>> pl = sx.plot_contour("time", [1], off_screen=True)

        The labels can be obtained using:

        >>> sx.get_all_label_spaces() # doctest: +ELLIPSIS
        [{'...': ..., '...': ...}, {'...': ..., '...': ...}]
        """
        self._evaluate_result()

        # check if complex label, not supported
        labels = self.result_fields_container.get_label_space(0)
        if DefinitionLabels.complex in labels.keys():
            raise core_errors.ComplexPlottingError

        # Get default values if only one field in the fields container
        if len(self.result_fields_container) == 1:
            lab_space = self.result_fields_container.get_label_space(0)
            display_option = [*lab_space][0]
            option_id = lab_space[display_option]

        # If plotting on a path
        if self._evaluator._path is not None:
            # Try and use the new DpfPlotter from PyDPF-Core
            try:
                from ansys.dpf.core.plotter import DpfPlotter as DpfPlotterObj
            except:
                raise dpf_errors.CoreVersionError(version='0.3.4')
            # Initialize the plotter
            pl = DpfPlotterObj(**kwargs)
            # Sort the fields according to options
            new_fields_container = self._sort_fields_container_with_labels(
                option_id, display_option)
            # Add each field with its associated mesh to the plotting
            for field_m in new_fields_container:
                mesh_m = field_m.meshed_region
                pl.add_field(field_m, mesh_m)
            # Add the mesh associated to the path
            pl.add_mesh(self._evaluator._model.metadata.meshed_region,
                        style="surface", show_edges=True, opacity=0.3, **kwargs)
            # Show the plot
            pl.show_figure(**kwargs)

        # If not plotting on a path
        else:
            # Initialize a Plotter
            pl = DpfPlotter(self._evaluator._model.metadata.meshed_region, **kwargs)
            # Create an equivalent field container
            if len(self.result_fields_container) == 1:
                fc = self.result_fields_container
            else:
                # sorts and creates a new fields_container with only the desired labels
                fc = self._sort_fields_container_with_labels(option_id, display_option)
            # Call Plotter.plot_contour (to change for use of DpfPlotter
            pl.plot_contour(fc, **kwargs)

    def _plot_chart(self):
        """Plot the minimum/maximum result values over time.

        This method works if the time_freq_support contains
        several time_steps (for example, a transient analysis).

        A time_scoping keyword must be used to select all the
        time_steps of the result.

        """
        self._evaluate_result()
        # tfq = self._evaluator._model.metadata.time_freq_support
        # timeids = list(range(1,tfq.n_sets+1))
        # res_op = self._evaluator._result_operator
        # res_op.inputs.time_scoping.connect(timeids)
        # new_fields_container = res_op.get_output(0, types.fields_container)
        pl = DpfPlotter(None)
        pl.plot_chart(self.result_fields_container)
