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

    @staticmethod
    def _get_component(component: str) -> str:
        """Returns a valid component for the operator naming."""
        component_input_map = {
            "X": "X",
            "Y": "Y",
            "Z": "Z",
            "0": "X",
            "1": "Y",
            "2": "Z",
            "XY": "XY",
            "YZ": "YZ",
            "XZ": "XZ",
            "3": "XY",
            "4": "YZ",
            "5": "XZ",
        }

        # Select the operator based on component
        comp = ""
        if component is not None:
            try:
                comp = component_input_map[component]
            except KeyError:
                ValueError(
                    f"Component {component} is not a valid component argument. "
                    f"Valid component arguments are:\n{list(component_input_map.keys())}"
                )
        return comp

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
                named_selection, server=self._model._server
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
        **kwargs,
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
        operator_name_base = "U"

        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False
        comp = self._get_component(component=str(component))
        op_name = operator_name_base + comp
        op = self._model.operator(name=op_name)
        wf.set_output_name("out", op.outputs.fields_container)
        wf.add_operator(op)

        time_scoping = self._select_time_freq(selection, steps)

        # Set the time_scoping if necessary
        if time_scoping:
            op.connect(0, time_scoping)

        # Build the mesh_scoping from nodes or selection
        mesh_scoping = self._select_mesh_scoping(
            selection, nodes, elements, named_selection
        )
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            op.connect(1, mesh_scoping)

        wf.add_operator(op)

        # We will use the DataObject thing here.
        wf.set_output_name("out", op.outputs.fields_container)

        columns = ["X", "Y", "Z"] if comp == "" else [comp]
        return DataObject(
            wf.get_output("out", core.types.fields_container),
            name=op_name,
            columns=columns,
            mesh_scoping=mesh_scoping,
        )

    def velocity(
        self,
        selection: Union[Selection, None] = None,
        steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        component: Union[int, str, List[str], None] = None,
        named_selection: Union[str, None] = None,
        **kwargs,
    ) -> DataObject:
        """Extract velocity results from the simulation.

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
        operator_name_base = "V"

        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False
        comp = self._get_component(component=str(component))
        op_name = operator_name_base + comp
        op = self._model.operator(name=op_name)
        wf.set_output_name("out", op.outputs.fields_container)
        wf.add_operator(op)

        time_scoping = self._select_time_freq(selection, steps)

        # Set the time_scoping if necessary
        if time_scoping:
            op.connect(0, time_scoping)

        # Build the mesh_scoping from nodes or selection
        mesh_scoping = self._select_mesh_scoping(
            selection, nodes, elements, named_selection
        )
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            op.connect(1, mesh_scoping)

        columns = ["X", "Y", "Z"] if comp == "" else [comp]
        return DataObject(
            wf.get_output("out", core.types.fields_container),
            name=op_name,
            columns=columns,
            mesh_scoping=mesh_scoping,
        )

    #     def nodal_stress(
    #         self,
    #         steps: Optional[list[int]] = None,
    #         components: Optional[Union[int, str, list[str]]] = None,
    #         nodes: Optional[list[int]] = None,
    #         elements: Optional[list[int]] = None,
    #         named_selection: Optional[str] = None,
    #         selection: Optional[Selection] = None,
    #         element_shape: Optional[core.elements._element_shapes] = None,
    #         ordered: bool = True,
    #         **kwargs
    #     ) -> ResultData:
    #         pass

    def elemental_stress(
        self,
        selection: Union[Selection, None] = None,
        steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        component: Union[int, str, List[str], None] = None,
        named_selection: Union[str, None] = None,
        element_shape: Union[
            core.elements._element_shapes, None
        ] = None,  # What is that for?
        **kwargs,
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

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        operator_name_base = "S"
        # TODO: What about Principal stresses?
        wf = core.Workflow(server=self._model._server)
        wf.progress_bar = False
        comp = self._get_component(component=str(component))
        op_name = operator_name_base + comp
        op = self._model.operator(name=op_name)
        wf.set_output_name("out", op.outputs.fields_container)
        wf.add_operator(op)

        time_scoping = self._select_time_freq(selection, steps)

        # Set the time_scoping if necessary
        if time_scoping:
            op.connect(0, time_scoping)

        # Build the mesh_scoping from nodes or selection
        mesh_scoping = self._select_mesh_scoping(
            selection, nodes, elements, named_selection
        )
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            op.connect(1, mesh_scoping)

        columns = ["X", "Y", "Z"] if comp == "" else [comp]
        return DataObject(
            wf.get_output("out", core.types.fields_container),
            name=op_name,
            columns=columns,
            mesh_scoping=mesh_scoping,
        )


#     def raw_stress(
#         self,
#         steps: Optional[list[int]] = None,
#         components: Optional[Union[int, str, list[str]]] = None,
#         nodes: Optional[list[int]] = None,
#         elements: Optional[list[int]] = None,
#         named_selection: Optional[str] = None,
#         selection: Optional[Selection] = None,
#         element_shape: Optional[core.elements._element_shapes] = None,
#         ordered: bool = True,
#         **kwargs
#     ) -> ResultData:
#         pass


class FluidSimulation(Simulation):
    """Provides a fluid type simulation."""

    def __init__(self, data_sources, model):
        """Initialize the class."""
        super().__init__(data_sources, model)

    def displacement(
        self,
        selection: Union[Selection, None] = None,  # Can deal with qualifiers
        steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        qualifiers: Union[List[str], None] = None,
        component: Union[int, str, List[str], None] = None,
        named_selection: Union[str, None] = None,
        **kwargs,
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
            qualifiers:
                List of qualifiers to get results for.
            component:
                Component to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        pass

    def velocity(
        self,
        selection: Union[Selection, None] = None,  # Can deal with qualifiers
        steps: Union[List[int], None] = None,
        nodes: Union[List[int], None] = None,
        elements: Union[List[int], None] = None,
        qualifiers: Union[List[str], None] = None,
        component: Union[int, str, List[str], None] = None,
        named_selection: Union[str, None] = None,
        **kwargs,
    ) -> DataObject:
        """Extract velocity results from the simulation.

        Args:
            selection:
                Selection to get results for.
            steps:
                List of steps to get results for.
            nodes:
                List of nodes to get results for.
            elements:
                List of elements to get results for.
            qualifiers:
                List of qualifiers to get results for.
            component:
                Component to get results for.
            named_selection:
                Named selection to get results for.

        Returns
        -------
            Returns a :class:`ansys.dpf.post.data_object.DataObject` instance.

        """
        pass
