from typing import Any, Union

from ansys.dpf.core import Workflow
from ansys.dpf.gate.common import locations

from ansys.dpf.post.misc import connect_any
from ansys.dpf.post.result_workflows._utils import _CreateOperatorCallable
from ansys.dpf.post.selection import _WfNames


def _create_averaging_workflow(
    has_skin: bool,
    location: str,
    mesh_averaging_needed: bool,
    create_operator_callable: _CreateOperatorCallable,
    server,
):
    average_wf = Workflow(server=server)

    input_data_fwd = create_operator_callable(name="forward_fc")
    averaged_data_fwd = create_operator_callable(name="forward_fc")

    mesh_averaging_input_fwd = create_operator_callable(name="forward_fc")
    average_wf.add_operators(
        [input_data_fwd, averaged_data_fwd, mesh_averaging_input_fwd]
    )
    mesh_averaging_input_fwd.connect(0, input_data_fwd, 0)

    average_wf.set_input_name(_WfNames.input_data, input_data_fwd)
    average_wf.set_output_name(_WfNames.output_data, averaged_data_fwd, 0)

    skin_mesh_fwd = create_operator_callable(name="forward")

    if has_skin:
        average_wf.add_operator(skin_mesh_fwd)
        average_wf.set_input_name(_WfNames.skin, skin_mesh_fwd)

        if server.meet_version("6.2"):
            solid_to_skin_operator = create_operator_callable(name="solid_to_skin_fc")
        else:
            solid_to_skin_operator = create_operator_callable(name="solid_to_skin")

        average_wf.add_operator(solid_to_skin_operator)
        solid_to_skin_operator.connect(0, input_data_fwd, 0)
        mesh_averaging_input_fwd.connect(0, solid_to_skin_operator, 0)

        if hasattr(solid_to_skin_operator.inputs, "mesh_scoping"):
            connect_any(solid_to_skin_operator.inputs.mesh_scoping, skin_mesh_fwd)
            # To keep for retro-compatibility
        else:
            connect_any(solid_to_skin_operator.inputs.mesh, skin_mesh_fwd)

        if server.meet_version("8.0"):
            # solid mesh_input only supported for server version
            # 8.0 and up
            average_wf.set_input_name(
                _WfNames.skin_input_mesh, solid_to_skin_operator.inputs.solid_mesh
            )

    if (
        location == locations.nodal
        or location == locations.elemental
        and mesh_averaging_needed
    ):
        if location == locations.nodal:
            operator_name = "to_nodal_fc"
        else:
            operator_name = "to_elemental_fc"
        mesh_average_op = create_operator_callable(name=operator_name)
        average_wf.add_operator(mesh_average_op)
        mesh_average_op.connect(0, mesh_averaging_input_fwd, 0)
        averaged_data_fwd.connect(0, mesh_average_op, 0)
    else:
        averaged_data_fwd.connect(0, mesh_averaging_input_fwd, 0)

    return average_wf


def _create_principal_workflow(
    components_to_extract: list[int],
    create_operator_callable: _CreateOperatorCallable,
    server,
):
    principal_wf = Workflow(server=server)

    # Instantiate the required operator
    principal_op = create_operator_callable(name="invariants_fc")
    principal_wf.add_operator(principal_op)
    principal_wf.set_input_name(_WfNames.input_data, principal_op)
    # Set as future output of the workflow
    if len(components_to_extract) == 1:
        principal_output = getattr(
            principal_op.outputs, f"fields_eig_{components_to_extract[0] + 1}"
        )
        principal_wf.set_output_name(_WfNames.output_data, principal_output)
    else:
        raise NotImplementedError("Cannot combine principal results yet.")
        # We need to define the behavior for storing different results in a DataFrame

    return principal_wf


def _create_equivalent_workflow(
    create_operator_callable: _CreateOperatorCallable, server
):
    equivalent_wf = Workflow(server=server)
    equivalent_op = create_operator_callable(name="eqv_fc")
    equivalent_wf.add_operator(operator=equivalent_op)
    equivalent_wf.set_input_name(_WfNames.input_data, equivalent_op)
    equivalent_wf.set_output_name(
        _WfNames.output_data, equivalent_op.outputs.fields_container
    )
    return equivalent_wf


def _create_mesh_workflow(mesh_provider: Any, server):
    mesh_wf = Workflow(server=server)
    mesh_wf.add_operator(mesh_provider)
    mesh_wf.set_output_name(_WfNames.initial_mesh, mesh_provider)
    return mesh_wf


def _create_extract_component_workflow(
    create_operator_callable: _CreateOperatorCallable,
    components_to_extract: list[int],
    base_name: str,
    server,
):
    extract_component_wf = Workflow(server=server)

    # Instantiate a component selector operator
    extract_op = create_operator_callable(name="component_selector_fc")
    # Feed it the current workflow output
    extract_component_wf.set_input_name(_WfNames.input_data, extract_op)

    # Feed it the requested components
    extract_op.connect(1, components_to_extract)
    extract_component_wf.add_operator(operator=extract_op)
    # Set as future output of the workflow
    extract_component_wf.set_output_name(
        _WfNames.output_data, extract_op.outputs.fields_container
    )

    result_is_single_component = False
    if len(components_to_extract) == 1:
        new_base_name = base_name + f"_{components_to_extract[0]}"
        result_is_single_component = True
    else:
        new_base_name = base_name

    return extract_component_wf, new_base_name, result_is_single_component


def _create_norm_workflow(
    create_operator_callable: _CreateOperatorCallable, base_name: str, server
):
    norm_wf = Workflow(server=server)
    norm_op = create_operator_callable(name="norm_fc")
    norm_wf.add_operator(operator=norm_op)
    norm_wf.set_input_name(_WfNames.input_data, norm_op)
    norm_wf.set_output_name(_WfNames.output_data, norm_op)
    new_base_name = base_name + "_N"
    return norm_wf, new_base_name


def _create_initial_result_workflow(
    name: str, server, create_operator_callable: _CreateOperatorCallable
):
    initial_result_workflow = Workflow(server=server)

    initial_result_op = create_operator_callable(name=name)
    initial_result_workflow.set_input_name(_WfNames.mesh, initial_result_op, 7)
    initial_result_workflow.set_input_name(_WfNames.location, initial_result_op, 9)

    initial_result_workflow.add_operator(initial_result_op)
    initial_result_workflow.set_output_name(_WfNames.output_data, initial_result_op, 0)
    initial_result_workflow.set_input_name(
        "time_scoping", initial_result_op.inputs.time_scoping
    )
    initial_result_workflow.set_input_name(
        "mesh_scoping", initial_result_op.inputs.mesh_scoping
    )

    initial_result_workflow.set_input_name(_WfNames.read_cyclic, initial_result_op, 14)
    initial_result_workflow.set_input_name(
        _WfNames.cyclic_sectors_to_expand, initial_result_op, 18
    )
    initial_result_workflow.set_input_name(_WfNames.cyclic_phase, initial_result_op, 19)

    return initial_result_workflow


def _create_sweeping_phase_workflow(
    create_operator_callable: _CreateOperatorCallable,
    server,
    amplitude: bool,
    sweeping_phase: Union[float, None],
):
    sweeping_phase_workflow = Workflow(server=server)
    # Add an optional sweeping phase or amplitude operation if requested
    # (must be after comp_selector for U)
    # (must be before norm operation for U)
    if sweeping_phase is not None and not amplitude:
        if isinstance(sweeping_phase, int):
            sweeping_phase = float(sweeping_phase)
        if not isinstance(sweeping_phase, float):
            raise ValueError("Argument sweeping_phase must be a float.")
        sweeping_op = create_operator_callable(name="sweeping_phase_fc")
        sweeping_op.connect(2, sweeping_phase)
        sweeping_op.connect(3, "degree")
        sweeping_op.connect(4, False)
        sweeping_phase_workflow.add_operator(operator=sweeping_op)

        sweeping_phase_workflow.set_input_name(_WfNames.input_data, sweeping_op)
        sweeping_phase_workflow.set_output_name(_WfNames.output_data, sweeping_op)
    elif amplitude:
        amplitude_op = create_operator_callable(name="amplitude_fc")
        sweeping_phase_workflow.add_operator(operator=amplitude_op)
        sweeping_phase_workflow.set_input_name(_WfNames.input_data, amplitude_op)
        sweeping_phase_workflow.set_output_name(_WfNames.output_data, amplitude_op)
    else:
        return None

    return sweeping_phase_workflow
