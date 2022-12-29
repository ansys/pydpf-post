"""Module containing the ``Simulation`` class."""
import re
from typing import List, Union

from ansys.dpf import core
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.selection import Selection


class Simulation:
    """Provides the main class of the DPF-Post solution."""

    def __init__(self, data_sources: core.DataSources, model: core.Model):
        """Initialize the solution using a ``dpf.core.Model`` object."""
        self._model = model
        self._data_sources = data_sources
        self._geometry = None
        self._active_selection = None
        self._mesh = None

    @property
    def results(self) -> List[str]:
        """Available results.

        Returns a list of available results as strings.
        """
        return self._model.metadata.result_info.available_results

    @property
    def mesh(self) -> Mesh:
        """Mesh representation of the model.

        Returns the :class:`ansys.dpf.post.mesh.Mesh` class.
        """
        if self._mesh is None:
            self._mesh = Mesh(self._model.metadata.meshed_region)
        return self._mesh

    def activate_selection(self, selection_object: Selection):
        """Selection currently active.

        Returns the current active :class:`ansys.dpf.post.selection.Selection` class.
        """
        self._active_selection = selection_object.selection

    def deactivate_selection(self):
        """Deactivate the currently active selection."""
        self._active_selection = None

    @property
    def _time_frequencies(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._model.metadata.time_freq_support

    @property
    def time_freq_support(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._time_frequencies

    def get_result_info(self):
        """Get result file information.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> solution = post.load_solution(examples.static_rst)
        >>> print(solution.get_result_info()) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        Static analysis
        Unit system: MKS: m, kg, N, s, V, A, degC
        Physics Type: ...
        Available results:
             -  displacement: Nodal Displacement
             -  reaction_force: Nodal Force
             -  stress: ElementalNodal Stress
             -  elemental_volume: Elemental Volume
             -  stiffness_matrix_energy: Elemental Energy-stiffness matrix
             -  artificial_hourglass_energy: Elemental Hourglass Energy
             -  thermal_dissipation_energy: Elemental thermal dissipation energy
             -  kinetic_energy: Elemental Kinetic Energy
             -  co_energy: Elemental co-energy
             -  incremental_energy: Elemental incremental energy
             -  elastic_strain: ElementalNodal Strain
             -  structural_temperature: ElementalNodal Temperature
        """
        return self._model.metadata.result_info

    def __str__(self):
        """Get the string representation of this class."""
        txt = (
            "%s." % re.sub(r"(?<!^)(?=[A-Z])", " ", type(self).__name__)
            + "\n\n\nData Sources\n------------------------------\n"
        )
        ds_str = self._model._data_sources.__str__()
        txt += ds_str
        txt += "\n\n"
        txt += self._model.__str__()
        return txt


class MechanicalSimulation(Simulation):
    """Provides a mechanical type solution."""

    def __init__(self, data_sources: core.DataSources, model: core.Model):
        """Instantiate a mechanical type solution."""
        self._columns = None
        super().__init__(data_sources, model)

    def _select_time_freq(self, selection=None, steps=None):
        """Select time."""
        # Build the time_scoping from steps or selection
        time_scoping = None
        if selection:
            time_scoping = selection.time_scoping  # needs to be changed.
        if steps:
            time_scoping = core.time_freq_scoping_factory.scoping_by_sets(
                steps, server=self._model._server
            )
        if not time_scoping:
            time_scoping = core.time_freq_scoping_factory.scoping_on_all_time_freqs(
                self._model.metadata.time_freq_support
            )
        return time_scoping

    def _select_mesh_scoping(
        self, selection=None, nodes=None, elements=None, named_selection=None
    ):
        if (nodes is not None or elements is not None) and named_selection is not None:
            raise ValueError(
                "nodes/elements and named_selection are mutually exclusive"
            )

        if selection is not None and (
            nodes is not None or named_selection is not None or elements is not None
        ):
            raise ValueError(
                "selection and nodes/elements/named_selection are mutually exclusive"
            )

        mesh_scoping = None
        if selection:
            mesh_scoping = selection.mesh_scoping

        if named_selection:
            mesh_scoping = core.mesh_scoping_factory.named_selection_scoping(
                named_selection, self._model, server=self._model._server
            )
        if nodes:
            mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
                nodes, server=self._model._server
            )

        if elements:
            mesh_scoping = core.mesh_scoping_factory.elemental_scoping(
                element_ids=elements, server=self._model._server
            )

        return mesh_scoping

    def displacement(
        self,
        selection: Union[Selection, None] = None,
        steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        component: Union[int, str, List[str], None] = None,
        named_selection: Union[str, None] = None,
        # ordered: bool = True,
        **kwargs
    ) -> DataObject:
        """Extract displacement results from the simulation.

        Args:
            selection:
                Selection to get results for.
            steps:
                List of steps to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            component:
                Component to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False

        # Select the operator based on component
        if component is None:
            op_name = "U"
            self._columns = ["X", "Y", "Z"]
        elif isinstance(component, (str, int)):
            if component in ["X", 0]:
                op_name = "UX"
                self._columns = ["X"]
            elif component in ["Y", 1]:
                op_name = "UY"
                self._columns = ["Y"]
            elif component in ["Z", 2]:
                op_name = "UZ"
                self._columns = ["Z"]
            else:
                op_name = "U"
                self._columns = ["X", "Y", "Z"]
        else:
            raise TypeError("Component must be a string or an integer")

        disp_op = self._model.operator(name=op_name)

        time_scoping = self._select_time_freq(selection, steps)

        # Set the time_scoping if specified. Return all times available if not.
        if time_scoping:
            disp_op.connect(0, time_scoping)

        # Build the mesh_scoping from nodes or selection
        mesh_scoping = self._select_mesh_scoping(
            selection, nodes, elements, named_selection
        )
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            disp_op.connect(1, mesh_scoping)

        wf.add_operator(disp_op)

        # We will use the DataObject thing here.
        wf.set_output_name("out", disp_op.outputs.fields_container)

        return DataObject(
            wf.get_output("out", core.types.fields_container),
            columns=self._columns,
            mesh_scoping=mesh_scoping,
        )

    # def velocity(
    #     self,
    #     selection: Union[Selection, None] = None,
    #     steps: Union[List[int], None] = None,
    #     nodes: Union[List[int], None] = None,
    #     elements: Union[List[int], None] = None,
    #     component: Union[int, str, List[str], None] = None,
    #     named_selection: Union[str, None] = None,
    #     # ordered: bool = True,
    #     **kwargs
    # ) -> DataObject:
    #     """Extract stress results from the simulation.

    #     Args:
    #         selection:
    #             Selection to get results for.
    #         steps:
    #             List of steps to get results for.
    #         nodes:
    #             List of nodes to get results for.
    #         elements:
    #             List of elements to get results for.
    #         component:
    #             Component to get results for.
    #         named_selection:
    #             Named selection to get results for.
    #         location:
    #             Stress location to get results for.

    #     Returns
    #     -------
    #         Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

    #     """
    #     wf = core.Workflow(server=self._model._server)
    #     wf.progress_bar = False

    #     # Select the operator based on component
    #     if component is None:
    #         op_name = "V"
    #     elif isinstance(component, (str, int)):
    #         if component in ["X", 0]:
    #             op_name = "VX"
    #         elif component in ["Y", 1]:
    #             op_name = "VY"
    #         elif component in ["Z", 2]:
    #             op_name = "VZ"
    #         else:
    #             op_name = "V"
    #     else:
    #         raise TypeError("Component must be a string or an integer")

    #     disp_op = self._model.operator(name=op_name)

    #     time_scoping = self._select_time_freq(selection, steps)

    #     # Set the time_scoping if necessary
    #     if time_scoping:
    #         disp_op.connect(0, time_scoping)

    #     # Build the mesh_scoping from nodes or selection
    #     mesh_scoping = self._select_mesh_scoping(
    #         selection, nodes, elements, named_selection
    #     )
    #     # Set the mesh_scoping if necessary
    #     if mesh_scoping:
    #         disp_op.connect(1, mesh_scoping)

    #     wf.add_operator(disp_op)

    #     # We will use the DataObject thing here.
    #     wf.set_output_name("out", disp_op.outputs.fields_container)

    #     return DataObject(
    #         wf.get_output("out", core.types.fields_container),
    #         mesh_scoping=mesh_scoping,
    #     )

    def nodal_stress(self, **kwargs):
        """Connect to the stress method with nodal location."""
        return self.stress(location="Nodal", **kwargs)

    def elemental_stress(self, **kwargs):
        """Connect to the stress method with elemental location."""
        return self.stress(location="Elemental", **kwargs)

    def raw_stress(self, **kwargs):
        """Connect to the stress method with elemental/nodal location."""
        return self.stress(location="ElementalNodal", **kwargs)

    def stress(
        self,
        location: str = "Nodal",
        selection: Union[Selection, None] = None,
        steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        component: Union[int, str, List[str], None] = None,
        named_selection: Union[str, None] = None,
        **kwargs
    ) -> DataObject:
        """Extract stress results from the simulation.

        Args:
            selection:
                Selection to get results for.
            steps:
                List of steps to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            component:
                Component to get results for.
            named_selection:
                Named selection to get results for.
            location:
                Stress location to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False

        # Select the operator based on component
        if component is None:
            op_name = "S"
            self._columns = ["X", "Y", "Z", "XY", "YZ", "XZ"]
        elif isinstance(component, str):
            if component == "X":
                op_name = "SX"
                self._columns = ["X"]
            elif component == "XY":
                op_name = "SXY"
                self._columns = ["XY"]
            elif component == "XZ":
                op_name = "SXZ"
                self._columns = ["XZ"]
            elif component == "Y":
                op_name = "SY"
                self._columns = ["Y"]
            elif component == "YZ":
                op_name = "SYZ"
                self._columns = ["YZ"]
            elif component == "Z":
                op_name = "SZ"
                self._columns = ["Z"]
            else:
                op_name = "S"
                self._columns = ["X", "Y", "Z", "XY", "YZ", "XZ"]
        else:
            raise TypeError("Component must be a string or an integer")

        disp_op = self._model.operator(name=op_name)
        disp_op.connect(9, location)

        time_scoping = self._select_time_freq(selection, steps)

        # Set the time_scoping if necessary
        if time_scoping:
            disp_op.connect(0, time_scoping)

        # Build the mesh_scoping from nodes or selection
        mesh_scoping = self._select_mesh_scoping(
            selection, nodes, elements, named_selection
        )
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            disp_op.connect(1, mesh_scoping)

        wf.add_operator(disp_op)

        # We will use the DataObject thing here.
        wf.set_output_name("out", disp_op.outputs.fields_container)

        return DataObject(
            wf.get_output("out", core.types.fields_container),
            columns=self._columns,
            mesh_scoping=mesh_scoping,
        )


# class FluidSolution(Simulation):
#     """Provides a fluid type solution."""

#     def __init__(self, data_sources, model):
#         super().__init__(data_sources, model)

#     def displacement(
#         self,
#         steps: Optional[list[int]] = None,
#         components: Optional[Union[int, str, list[str]]] = None,
#         nodes: Optional[list[int]] = None,
#         named_selection: Optional[str] = None,
#         selection: Optional[Selection] = None,
#         ordered: bool = True,
#         **kwargs
#     ) -> ResultData:
#         pass

#     def velocity(
#         self,
#         steps: Optional[list[int]] = None,
#         components: Optional[Union[int, str, list[str]]] = None,
#         nodes: Optional[list[int]] = None,
#         named_selection: Optional[str] = None,
#         selection: Optional[Selection] = None,  # Can deal with qualifiers
#         ordered: bool = True,
#         **kwargs
#     ) -> ResultData:
#         pass
