"""Module containing the ``Simulation`` class."""
from abc import ABC, abstractmethod
from enum import Enum
import re
from typing import List, Union

from ansys.dpf.core import DataSources, Model
from ansys.dpf.core.plotter import DpfPlotter

from ansys.dpf import core
from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.selection import Selection


class ResultCategory(Enum):
    """Enum for available result categories."""

    scalar = 1
    vector = 2
    matrix = 3
    principal = 4
    equivalent = 5


class Simulation(ABC):
    """Base class of all PyDPF-Post simulation types."""

    _component_to_id = {
        "1": 0,
        "2": 1,
        "3": 2,
        "4": 3,
        "5": 4,
        "6": 5,
        "X": 0,
        "Y": 1,
        "Z": 2,
        "XY": 3,
        "YZ": 4,
        "XZ": 5,
    }

    _component_names = ["X", "Y", "Z", "XX", "XY", "YZ"]
    _principal_names = ["1", "2", "3"]

    def __init__(self, data_sources: DataSources, model: Model):
        """Initialize the simulation using a ``dpf.core.Model`` object."""
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
        r"""Available results.

        Returns a list of available results as strings.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.results) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['displacement\nOperator name: "U"\n...Units: degc\n']
        """
        return [
            str(result) for result in self._model.metadata.result_info.available_results
        ]

    @property
    def geometries(self):
        """List of constructed geometries in the simulation.

        Returns a list of geometry objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.geometries) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._geometries

    @property
    def boundary_conditions(self):
        """List of boundary conditions in the simulation.

        Returns a list of boundary_condition objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.boundary_conditions) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._boundary_conditions

    @property
    def loads(self):
        """List of loads in the simulation.

        Returns a list of load objects.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.loads) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        return self._loads

    @property
    def mesh(self) -> Mesh:
        """Mesh representation of the model.

        Returns a :class:`ansys.dpf.post.mesh.Mesh` object.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.mesh) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.mesh.Mesh object at ...>
        """
        if self._mesh is None:
            self._mesh = Mesh(self._model.metadata.meshed_region)
        return self._mesh

    @property
    def named_selections(self) -> List[str]:
        """List of named selections in the simulation.

        Returns a list of named selections names.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> print(simulation.named_selections) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['_FIXEDSU']
        """
        if self._named_selections is None:
            self._named_selections = self._model.metadata.available_named_selections
        return self._named_selections

    def plot(
        self,
        mesh: bool = True,
        constructed_geometries: bool = True,
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
            constructed_geometries:
                Whether to plot the constructed geometries.
            loads:
                Whether to plot the loads.
            boundary_conditions:
                Whether to plot the boundary conditions.

        Returns
        -------
            Returns a plotter instance of the active visualization backend.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> simulation.plot() # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        """
        plt = DpfPlotter()
        if mesh:
            plt.add_mesh(self.mesh._meshed_region)
        if constructed_geometries:
            for geom in self.geometries:
                getattr(plt, "add_" + str(type(geom).__name__).lower())(geom)
        if loads:
            pass
        if boundary_conditions:
            pass
        plt.show_figure()

    @property
    def active_selection(self) -> Selection:
        """Active selection used by default for result queries.

        Returns a :object:`ansys.dpf.post.selection.Selection` object.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.activate_selection(selection=selection)
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        """
        return self._active_selection

    def activate_selection(self, selection: Selection):
        """Sets a selection as active on the simulation.

        Activating a given selection on a simulation means it is used
        as a default selection/filter in further result queries.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.activate_selection(selection=selection)
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        """
        self._active_selection = selection

    def deactivate_selection(self):
        """Deactivate the currently active selection.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> simulation = post.load_simulation(examples.static_rst)
        >>> selection = post.selection.Selection()
        >>> simulation.activate_selection(selection=selection)
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        <ansys.dpf.post.selection.Selection object at ...>
        >>> simulation.deactivate_selection()
        >>> print(simulation.active_selection) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        None
        """
        self._active_selection = None

    @property
    def _time_frequencies(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._model.metadata.time_freq_support

    @property
    def time_freq_support(self):
        """Description of the temporal/frequency analysis of the model."""
        return self._time_frequencies

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

    def _build_components_from_components(self, base_name, category, components):
        # Create operator internal names based on components
        out = []
        if components is None:
            out = None
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if not isinstance(components, list):
                raise ValueError(
                    "Argument 'components' must be an int, a str, or a list of either."
                )
            for comp in components:
                if not (isinstance(comp, str) or isinstance(comp, int)):
                    raise ValueError(
                        "Argument 'components' can only contain integers and/or strings."
                    )
                if isinstance(comp, int):
                    comp = str(comp)
                if comp not in self._component_to_id.keys():
                    raise ValueError(
                        f"Component {comp} is not valid. Please use one of: "
                        f"{list(self._component_to_id.keys())}."
                    )
                out.append(self._component_to_id[comp])

        # Take unique values and build names list
        if out is not None:
            out = list(set(out))
        if out is None and category == ResultCategory.vector:
            columns = [base_name + comp for comp in self._component_names[:3]]
        elif out is None and category == ResultCategory.matrix:
            columns = [base_name + comp for comp in self._component_names]
        else:
            columns = [base_name + self._component_names[i] for i in out]
        return out, columns

    def _build_components_from_principal(self, base_name, components):
        # Create operator internal names based on principal components
        out = []
        if components is None:
            out = None
        else:
            if isinstance(components, int) or isinstance(components, str):
                components = [components]
            if not isinstance(components, list):
                raise ValueError(
                    "Argument 'components' must be an int, a str, or a list of either."
                )
            for comp in components:
                if not (isinstance(comp, str) or isinstance(comp, int)):
                    raise ValueError(
                        "Argument 'components' can only contain integers and/or strings."
                    )
                if str(comp) not in self._principal_names:
                    raise ValueError(
                        "A principal component ID must be one of: "
                        f"{self._principal_names}."
                    )
                out.append(comp - 1)

        # Take unique values
        if out is not None:
            out = list(set(out))
        # Build columns names
        if out is None:
            columns = [base_name + str(comp) for comp in self._principal_names]
        else:
            columns = [base_name + self._principal_names[i] for i in out]
        return out, columns

    @abstractmethod
    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass

    @abstractmethod
    def _build_mesh_scoping(self) -> core.mesh_scoping_factory.Scoping:
        """Generate a mesh_scoping from input arguments."""
        pass

    def _build_result_operator(self, name, time_scoping, mesh_scoping, location):
        op = self._model.operator(name=name)
        # Set the time_scoping if necessary
        if time_scoping:
            op.connect(0, time_scoping)
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            op.connect(1, mesh_scoping)

        op.connect(7, self.mesh._meshed_region)
        op.connect(9, location)
        return op


class MechanicalSimulation(Simulation, ABC):
    """Base class for mechanical type simulations.

    This class provides common methods and properties for all mechanical type simulations.
    """

    def __init__(self, data_sources: core.DataSources, model: core.Model):
        """Instantiate a mechanical type simulation."""
        super().__init__(data_sources, model)

    def _build_mesh_scoping(
        self,
        selection=None,
        nodes=None,
        elements=None,
        named_selections=None,
        location=core.locations.nodal,
    ) -> Union[core.mesh_scoping_factory.Scoping, None]:
        if (nodes is not None or elements is not None) and named_selections is not None:
            raise ValueError(
                "nodes/elements and named_selection are mutually exclusive"
            )

        if selection is not None and (
            nodes is not None or named_selections is not None or elements is not None
        ):
            raise ValueError(
                "selection and nodes/elements/named_selection are mutually exclusive"
            )

        mesh_scoping = None
        if selection:
            mesh_scoping = selection.mesh_scoping

        if named_selections:
            if type(named_selections) == str:
                mesh_scoping = core.mesh_scoping_factory.named_selection_scoping(
                    named_selections, server=self._model._server, model=self._model
                )
            elif type(named_selections) == list:
                merge_scopings_op = self._model.operator(name="merge::scoping")
                for pin, named_selection in enumerate(named_selections):
                    mesh_scoping_on_ns_op = self._model.operator(
                        name="scoping_provider_by_ns"
                    )
                    mesh_scoping_on_ns_op.connect(0, location)
                    mesh_scoping_on_ns_op.connect(1, named_selection)
                    merge_scopings_op.connect(
                        pin, mesh_scoping_on_ns_op.outputs.mesh_scoping
                    )
                mesh_scoping = merge_scopings_op.eval()

        if nodes:
            mesh_scoping = core.mesh_scoping_factory.nodal_scoping(
                nodes, server=self._model._server
            )

        if elements:
            mesh_scoping = core.mesh_scoping_factory.elemental_scoping(
                element_ids=elements, server=self._model._server
            )

        if mesh_scoping is None:
            return mesh_scoping

        # Transpose location if necessary
        if (
            location == core.locations.nodal
            and mesh_scoping.location != core.locations.nodal
        ) or (
            location == core.locations.elemental
            and mesh_scoping.location != core.locations.elemental
        ):
            mesh_scoping = core.operators.scoping.transpose(
                mesh_scoping=mesh_scoping,
                meshed_region=self.mesh._meshed_region,
                inclusive=1,
            ).eval()

        return mesh_scoping


class ModalMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical modal simulations."""

    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass


class HarmonicMechanicalSimulation(MechanicalSimulation):
    """Provides methods for mechanical harmonic simulations."""

    def _build_time_freq_scoping(self) -> core.time_freq_scoping_factory.Scoping:
        """Generate a time_freq_scoping from input arguments."""
        pass
