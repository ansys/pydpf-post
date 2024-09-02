import dataclasses
from typing import Any, Callable, List, Optional, Union

from ansys.dpf.core import Operator, Workflow
from ansys.dpf.core.available_result import _result_properties
from ansys.dpf.gate.common import locations

from ansys.dpf.post.result_workflows._component_helper import (
    ResultCategory,
    _create_components,
)
from ansys.dpf.post.result_workflows._sub_workflows import (
    _create_averaging_workflow,
    _create_equivalent_workflow,
    _create_extract_component_workflow,
    _create_initial_result_workflow,
    _create_mesh_workflow,
    _create_norm_workflow,
    _create_principal_workflow,
    _create_sweeping_phase_workflow,
)
from ansys.dpf.post.result_workflows._utils import _CreateOperatorCallable
from ansys.dpf.post.selection import Selection, _WfNames


@dataclasses.dataclass
class ResultWorkflows:
    """Contains all the sub-workflows needed to compute a result and some additional information.

    Use _create_result_workflows to create this object.
    Some workflows are optional. If they are not needed for a particular result
    , they are set to None.
    """

    # Workflow to compute the initial result (e.g. stress, displacement, etc.)
    # The location of this result is always elemental_nodal if force_elemental_nodal is True
    initial_result_workflow: Workflow
    # Workflow to average the result. Maps results on skin if needed and averages to the requested
    # location if force_elemental_nodal is True
    averaging_workflow: Workflow
    # The name of the requested result operator with
    # some modifications (e.g. "_VM" for equivalent stress)
    base_name: str
    # If True the initial_result_workflow requests the result at the elemental nodal location
    # and the averaging_workflow averages the result to the requested location. This is the
    # case for instance for skin results.
    force_elemental_nodal: bool
    # If True, the equivalent_workflow is computed before the averaging_workflow
    compute_equivalent_before_average: bool = False
    # List of component names at the end of the workflow. If None, the result is a scalar.
    components: Optional[list[str]] = None
    # Workflow to compute the mesh
    mesh_workflow: Optional[Workflow] = None
    # Workflow to compute the principal components of the result
    principal_workflow: Optional[Workflow] = None
    # Workflow to compute the equivalent result
    equivalent_workflow: Optional[Workflow] = None
    # Workflow normalize the result
    norm_workflow: Optional[Workflow] = None
    # Workflow to extract components of the result
    component_extraction_workflow: Optional[Workflow] = None
    # Workflow to sweep the phase of the result
    sweeping_phase_workflow: Optional[Workflow] = None


@dataclasses.dataclass
class _AveragingWorkflowInputs:
    location: Union[locations, str]
    force_elemental_nodal: bool


@dataclasses.dataclass
class _MeshWorkflowInputs:
    mesh_provider: Any


@dataclasses.dataclass
class _SweepingPhaseWorkflowInputs:
    amplitude: bool = (False,)
    sweeping_phase: Union[float, None] = (None,)


@dataclasses.dataclass
class _CreateWorkflowInputs:
    averaging_workflow_inputs: _AveragingWorkflowInputs
    has_skin: bool
    has_equivalent: bool
    has_principal: bool
    has_norm: bool
    base_name: str
    component_names: list[str]
    components_to_extract: list[int]
    should_extract_components: bool
    mesh_workflow_inputs: Optional[_MeshWorkflowInputs] = None
    sweeping_phase_workflow_inputs: Optional[_SweepingPhaseWorkflowInputs] = None


def _requires_manual_averaging(
    base_name: str,
    location: str,
    category: ResultCategory,
    selection: Optional[Selection],
    create_operator_callable: Callable[[str], Operator],
):
    res = _result_properties[base_name] if base_name in _result_properties else None
    if category == ResultCategory.equivalent and base_name[0] == "E":  # strain eqv
        return True
    if res is not None and selection is not None:
        return selection.requires_manual_averaging(
            location=location,
            result_native_location=res["location"],
            is_model_cyclic=create_operator_callable("is_cyclic").eval(),
        )
    return False


def _create_result_workflows(
    server,
    create_operator_callable: _CreateOperatorCallable,
    create_workflow_inputs: _CreateWorkflowInputs,
) -> ResultWorkflows:
    """Internal function to create the result workflows. Use _create_result_workflows instead."""
    initial_result_wf = _create_initial_result_workflow(
        name=create_workflow_inputs.base_name,
        server=server,
        create_operator_callable=create_operator_callable,
    )

    force_elemental_nodal = (
        create_workflow_inputs.averaging_workflow_inputs.force_elemental_nodal
    )
    average_wf = _create_averaging_workflow(
        location=create_workflow_inputs.averaging_workflow_inputs.location,
        has_skin=create_workflow_inputs.has_skin,
        force_elemental_nodal=force_elemental_nodal,
        create_operator_callable=create_operator_callable,
        server=server,
    )

    result_workflows: ResultWorkflows = ResultWorkflows(
        initial_result_workflow=initial_result_wf,
        averaging_workflow=average_wf,
        base_name=create_workflow_inputs.base_name,
        force_elemental_nodal=force_elemental_nodal,
        components=create_workflow_inputs.component_names,
    )

    if create_workflow_inputs.mesh_workflow_inputs is not None:
        result_workflows.mesh_workflow = _create_mesh_workflow(
            mesh_provider=create_workflow_inputs.mesh_workflow_inputs.mesh_provider,
            server=server,
        )

    if create_workflow_inputs.has_principal:
        result_workflows.principal_workflow = _create_principal_workflow(
            components_to_extract=create_workflow_inputs.components_to_extract,
            create_operator_callable=create_operator_callable,
            server=server,
        )

    if create_workflow_inputs.has_equivalent:
        result_workflows.equivalent_workflow = _create_equivalent_workflow(
            create_operator_callable=create_operator_callable, server=server
        )
        result_workflows.base_name += "_VM"
        # equivalent computation is done before averaging for strain because Mechanical
        # does it this way (MAPDL has a result named EPEL_EQV in the rst which
        # Mechanical uses directly
        if create_workflow_inputs.base_name[0] == "E":
            result_workflows.compute_equivalent_before_average = True

    if create_workflow_inputs.should_extract_components:
        (
            extract_component_wf,
            base_name,
            result_is_single_component,
        ) = _create_extract_component_workflow(
            create_operator_callable=create_operator_callable,
            components_to_extract=create_workflow_inputs.components_to_extract,
            component_names=create_workflow_inputs.component_names,
            base_name=create_workflow_inputs.base_name,
            server=server,
        )
        result_workflows.component_extraction_workflow = extract_component_wf
        if result_is_single_component:
            result_workflows.components = None
        result_workflows.base_name = base_name

    if create_workflow_inputs.has_norm:
        norm_wf, base_name = _create_norm_workflow(
            create_operator_callable=create_operator_callable,
            base_name=create_workflow_inputs.base_name,
            server=server,
        )
        result_workflows.norm_workflow = norm_wf
        result_workflows.components = None
        result_workflows.base_name = base_name

    if create_workflow_inputs.sweeping_phase_workflow_inputs is not None:
        result_workflows.sweeping_phase_workflow = _create_sweeping_phase_workflow(
            create_operator_callable=create_operator_callable,
            server=server,
            amplitude=create_workflow_inputs.sweeping_phase_workflow_inputs.amplitude,
            sweeping_phase=create_workflow_inputs.sweeping_phase_workflow_inputs.sweeping_phase,
        )

    return result_workflows


def _create_result_workflow_inputs(
    base_name: str,
    category: ResultCategory,
    components: Union[str, List[str], int, List[int], None],
    location: str,
    norm: bool,
    selection: Selection,
    create_operator_callable: Callable[[str], Operator],
    mesh_provider: Any,
    amplitude: bool = False,
    sweeping_phase: Union[float, None] = 0.0,
) -> _CreateWorkflowInputs:
    """Creates a CreateWorkflowInputs object to be used to create the result workflows."""
    component_names, components_to_extract, _ = _create_components(
        base_name, category, components
    )

    force_elemental_nodal = _requires_manual_averaging(
        base_name=base_name,
        location=location,
        category=category,
        selection=selection,
        create_operator_callable=create_operator_callable,
    )

    averaging_workflow_inputs = _AveragingWorkflowInputs(
        location=location,
        force_elemental_nodal=force_elemental_nodal,
    )

    mesh_workflow_inputs: Optional[_MeshWorkflowInputs] = None
    if selection.spatial_selection.requires_mesh:
        mesh_workflow_inputs = _MeshWorkflowInputs(mesh_provider=mesh_provider)

    has_principal = category == ResultCategory.principal

    should_extract_components = (
        category in [ResultCategory.vector, ResultCategory.matrix]
    ) and components_to_extract is not None

    sweeping_phase_workflow_inputs: Optional[_SweepingPhaseWorkflowInputs] = None
    if amplitude or sweeping_phase is not None:
        sweeping_phase_workflow_inputs = _SweepingPhaseWorkflowInputs(
            amplitude=amplitude,
            sweeping_phase=sweeping_phase,
        )

    return _CreateWorkflowInputs(
        base_name=base_name,
        averaging_workflow_inputs=averaging_workflow_inputs,
        has_skin=_WfNames.skin in selection.spatial_selection._selection.output_names,
        has_norm=norm,
        component_names=component_names,
        components_to_extract=components_to_extract,
        should_extract_components=should_extract_components,
        has_principal=has_principal,
        mesh_workflow_inputs=mesh_workflow_inputs,
        has_equivalent=category == ResultCategory.equivalent,
        sweeping_phase_workflow_inputs=sweeping_phase_workflow_inputs,
    )
