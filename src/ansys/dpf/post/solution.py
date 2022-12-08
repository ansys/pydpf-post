"""Module containing the ``Solution`` class.
"""
import re
from typing import Union, Optional

from ansys.dpf.post.mesh import Mesh
from ansys.dpf.post.selection import Selection
from ansys.dpf.post.result_data import ResultData

from ansys.dpf import core


class Solution:
    """Provides the main class of the DPF-Post solution."""

    def __init__(self, model: core.Model):
        """Initialize the solution using a ``dpf.core.Model`` object."""
        self._model = model
        self._active_selection = None
        self._mesh = None

    @property
    def mesh(self) -> Mesh:
        """Mesh representation of the model.

        Returns the :class:`ansys.dpf.post.mesh.Mesh` class.
        """
        if self._mesh is None:
            self._mesh = Mesh(self._model.metadata.meshed_region)
        return self._mesh

    def activate_selection(self, selection_object: Selection):
        self._active_selection = selection_object.selection

    def deactivate_selection(self):
        self._active_selection = None

    @property
    def _time_frequencies(self):
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


class MechanicalSolution(Solution):
    """Provides a mechanical type solution."""

    def __init__(self, model: core.Model):
        super().__init__(model)

    def displacement(self,
                     steps: Optional[list[int]] = None,
                     component: Optional[Union[int, str, list[str]]] = None,
                     nodes: Optional[list[int]] = None,
                     named_selection: Optional[str] = None,
                     selection: Optional[Selection] = None,
                     ordered: bool = True,
                     **kwargs
                     ) -> ResultData:
        wf = core.Workflow(server=self._model._server)

        # Select the operator based on component
        if component in ['X', 0]:
            op_name = "UX"
        elif component in ['Y', 1]:
            op_name = "UY"
        elif component in ['Z', 2]:
            op_name = "UZ"
        else:
            op_name = "U"
        disp_op = self._model.operator(name=op_name)

        # Build the time_scoping from steps or selection
        time_scoping = None
        if selection:
            time_scoping = selection.time_scoping
        if steps:
            time_scoping = core.time_freq_scoping_factory.scoping_by_sets(steps,
                                                                          server=self._model._server)
        # Set the time_scoping if necessary
        if time_scoping:
            disp_op.connect(0, time_scoping)

        # Build the mesh_scoping from nodes or selection
        mesh_scoping = None
        if selection:
            mesh_scoping = selection.mesh_scoping
        if named_selection:
            mesh_scoping = core.mesh_scoping_factory.named_selection_scoping(named_selection,
                                                                             server=self._model._server)
        if nodes:
            mesh_scoping = core.mesh_scoping_factory.nodal_scoping(nodes,
                                                                   server=self._model._server)
        # Set the mesh_scoping if necessary
        if mesh_scoping:
            disp_op.connect(1, mesh_scoping)

        wf.add_operator(disp_op)

        if ordered:
            ord_op = self._model.operator(name="Rescope_fc")
            ord_op.connect(0, disp_op.get_output(0))
            wf.set_output_name("out", ord_op.get_output(0))
        else:
            wf.set_output_name("out", disp_op.get_output(0))

        fc = wf.get_output("out", core.types.fields_container)
        return ResultData(fc)

    def velocity(self,
                 steps: Optional[list[int]] = None,
                 components: Optional[Union[int, str, list[str]]] = None,
                 nodes: Optional[list[int]] = None,
                 named_selection: Optional[str] = None,
                 selection: Optional[Selection] = None,
                 ordered: bool = True,
                 **kwargs
                 ) -> ResultData:
        pass

    def nodal_stress(self,
                     steps: Optional[list[int]] = None,
                     components: Optional[Union[int, str, list[str]]] = None,
                     nodes: Optional[list[int]] = None,
                     elements: Optional[list[int]] = None,
                     named_selection: Optional[str] = None,
                     selection: Optional[Selection] = None,
                     element_shape: Optional[core.elements._element_shapes] = None,
                     ordered: bool = True,
                     **kwargs
                     ) -> ResultData:
        pass

    def elemental_stress(self,
                         steps: Optional[list[int]] = None,
                         components: Optional[Union[int, str, list[str]]] = None,
                         nodes: Optional[list[int]] = None,
                         elements: Optional[list[int]] = None,
                         named_selection: Optional[str] = None,
                         selection: Optional[Selection] = None,
                         element_shape: Optional[core.elements._element_shapes] = None,
                         ordered: bool = True,
                         **kwargs
                         ) -> ResultData:
        pass

    def raw_stress(self,
                   steps: Optional[list[int]] = None,
                   components: Optional[Union[int, str, list[str]]] = None,
                   nodes: Optional[list[int]] = None,
                   elements: Optional[list[int]] = None,
                   named_selection: Optional[str] = None,
                   selection: Optional[Selection] = None,
                   element_shape: Optional[core.elements._element_shapes] = None,
                   ordered: bool = True,
                   **kwargs
                   ) -> ResultData:
        pass


class FluidSolution(Solution):
    """Provides a fluid type solution."""

    def __init__(self, data_sources, model):
        super().__init__(data_sources, model)

    def displacement(self,
                     steps: Optional[list[int]] = None,
                     components: Optional[Union[int, str, list[str]]] = None,
                     nodes: Optional[list[int]] = None,
                     named_selection: Optional[str] = None,
                     selection: Optional[Selection] = None,
                     ordered: bool = True,
                     **kwargs
                     ) -> ResultData:
        pass

    def velocity(self,
                 steps: Optional[list[int]] = None,
                 components: Optional[Union[int, str, list[str]]] = None,
                 nodes: Optional[list[int]] = None,
                 named_selection: Optional[str] = None,
                 selection: Optional[Selection] = None,  # Can deal with qualifiers
                 ordered: bool = True,
                 **kwargs
                 ) -> ResultData:
        pass
