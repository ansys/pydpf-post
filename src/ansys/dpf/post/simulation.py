"""Module containing the ``Simulation`` class."""
import re
from typing import List, Union

from ansys.dpf.core import DataSources, Model
from ansys.dpf.core.plotter import DpfPlotter

from ansys.dpf import core
from ansys.dpf.post.data_object import DataObject
from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.selection import Selection


class Simulation:
    """Provides the main class of the DPF-Post solution."""

    def __init__(self, data_sources: DataSources, model: Model):
        """Initialize the solution using a ``dpf.core.Model`` object."""
        self._model = model
        self._data_sources = data_sources
        self._geometries = []
        self._boundary_conditions = []
        self._loads = []
        self._active_selection = None
        self._named_selections = None
        self._mesh = None

    @property
    def results(self) -> List[str]:
        """Available results.

        Returns a list of available results as strings.
        """
        return [
            str(result) for result in self._model.metadata.result_info.available_results
        ]

    @property
    def geometries(self):
        """List of constructed geometries in the simulation.

        Returns a list of geometry objects.
        """
        return self._geometries

    @property
    def boundary_conditions(self):
        """List of boundary conditions in the simulation.

        Returns a list of boundary_condition objects.
        """
        return self._boundary_conditions

    @property
    def loads(self):
        """List of loads in the simulation.

        Returns a list of load objects.
        """
        return self._loads

    @property
    def mesh(self) -> Mesh:
        """Mesh representation of the model.

        Returns a :class:`ansys.dpf.post.mesh.Mesh` object.
        """
        if self._mesh is None:
            self._mesh = Mesh(self._model.metadata.meshed_region)
        return self._mesh

    @property
    def named_selections(self) -> List[str]:
        """List of named selections in the simulation.

        Returns a list of named selections names.
        """
        if self._named_selections is None:
            self._named_selections = self._model.metadata.available_named_selections
        return self._named_selections

    def plot(
        self,
        mesh: bool = True,
        geometry: bool = True,
        loads: bool = True,
        boundary_conditions: bool = True,
    ):
        """General plot of the simulation object.

        Plots by default the complete mesh contained in the simulation,
        as well as a representation of the constructed geometry,
        the loads, and the boundary conditions currently defined.
        Each representation can be deactivated with its respective boolean argument.

        Args:
            mesh:
                Whether to plot the mesh representation.
            geometry:
                Whether to plot the geometries.
            loads:
                Whether to plot the loads.
            boundary_conditions:
                Whether to plot the boundary conditions.

        Returns
        -------
            Returns a plotter instance of the active visualization backend.
        """
        plt = DpfPlotter()
        if mesh:
            plt.add_mesh(self.mesh._meshed_region)
        if geometry:
            for geom in self.geometries:
                getattr(plt, "add_" + str(type(geom).__name__).lower())(geom)
        if loads:
            pass
        if boundary_conditions:
            pass
        plt.show_figure()

    def activate_selection(self, selection: Selection):
        """Sets a selection as active on the simulation.

        Activating a given selection on a simulation means it is used
        as a default selection/filter in further result queries.
        """
        self._active_selection = selection

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
        # ordered: bool = True,
        **kwargs
    ) -> DataObject:
        """Extract displacement results from the solution.

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
        elif isinstance(component, (str, int)):
            if component in ["X", 0]:
                op_name = "UX"
            elif component in ["Y", 1]:
                op_name = "UY"
            elif component in ["Z", 2]:
                op_name = "UZ"
            else:
                op_name = "U"
        else:
            raise TypeError("Component must be a string or an integer")

        disp_op = self._model.operator(name=op_name)

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

        # Reorder
        # ord_op = self._model.operator(name="Rescope_fc")
        # ord_op.inputs.fields_container.connect(disp_op.outputs.fields_container)
        # ord_op.inputs.mesh_scoping.connect(mesh_scoping)

        # ord_op.connect(0, disp_op.outputs.fields_container)

        # We will use the DataObject thing here.
        wf.set_output_name("out", disp_op.outputs.fields_container)

        return DataObject(
            wf.get_output("out", core.types.fields_container),
            columns=["X", "Y", "Z"],
            mesh_scoping=mesh_scoping,
        )


#     def velocity(
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

#     def elemental_stress(
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


class StaticMechanicalSimulation(MechanicalSimulation):
    """Provides a mechanical static type simulation."""

    @property
    def steps(self) -> List[float]:
        """List of the steps in the simulation."""
        return self._time_frequencies.time_frequencies.data


class TransientMechanicalSimulation(MechanicalSimulation):
    """Provides a mechanical transient type simulation."""


class ModalMechanicalSimulation(MechanicalSimulation):
    """Provides a mechanical modal type simulation."""


class HarmonicMechanicalSimulation(MechanicalSimulation):
    """Provides a mechanical harmonic type simulation."""


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
