import dataclasses
from typing import Any, Callable, Optional, Protocol, Union

from ansys.dpf.core import MeshedRegion, Operator, Scoping, ScopingsContainer, Workflow
from ansys.dpf.gate.common import locations

from ansys.dpf.post.misc import connect_any
from ansys.dpf.post.selection import Selection, _WfNames


class CreateOperatorCallable(Protocol):
    def __call__(self, name: str) -> Operator:
        ...


def create_averaging_workflow(
    has_skin: bool,
    location: str,
    mesh_averaging_needed: bool,
    create_operator_callable: CreateOperatorCallable,
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


def create_principal_workflow(
    components_to_extract: list[int],
    create_operator_callable: CreateOperatorCallable,
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


def create_equivalent_workflow(
    create_operator_callable: CreateOperatorCallable, server
):
    equivalent_wf = Workflow(server=server)
    equivalent_op = create_operator_callable(name="eqv_fc")
    equivalent_wf.add_operator(operator=equivalent_op)
    equivalent_wf.set_input_name(_WfNames.input_data, equivalent_op)
    equivalent_wf.set_output_name(
        _WfNames.output_data, equivalent_op.outputs.fields_container
    )
    return equivalent_wf


def create_mesh_workflow(mesh_provider: Any, server):
    mesh_wf = Workflow(server=server)
    mesh_wf.add_operator(mesh_provider)
    mesh_wf.set_output_name(_WfNames.initial_mesh, mesh_provider)
    return mesh_wf


def create_extract_component_workflow(
    create_operator_callable: CreateOperatorCallable,
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


def create_norm_workflow(
    create_operator_callable: CreateOperatorCallable, base_name: str, server
):
    norm_wf = Workflow(server=server)
    norm_op = create_operator_callable(name="norm_fc")
    norm_wf.add_operator(operator=norm_op)
    norm_wf.set_input_name(_WfNames.input_data, norm_op)
    norm_wf.set_output_name(_WfNames.output_data, norm_op)
    new_base_name = base_name + "_N"
    return norm_wf, new_base_name


def create_output_workflow(create_operator_callable: CreateOperatorCallable, server):
    forward_output_wf = Workflow(server=server)
    forward_op = create_operator_callable(name="forward")
    forward_output_wf.add_operator(forward_op)
    forward_output_wf.set_input_name(_WfNames.input_data, forward_op)
    forward_output_wf.set_output_name("out", forward_op)
    return forward_output_wf


def create_initial_result_workflow(
    name: str,
    location: Union[locations, str],
    force_elemental_nodal: bool,
    server,
    create_operator_callable: CreateOperatorCallable,
    mesh: MeshedRegion,
):
    initial_result_workflow = Workflow(server=server)

    initial_result_op = create_operator_callable(name=name)
    initial_result_op.connect(7, mesh)
    if force_elemental_nodal:
        initial_result_op.connect(9, "ElementalNodal")
    elif location:
        initial_result_op.connect(9, location)

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


def create_sweeping_phase_workflow(
    create_operator_callable: CreateOperatorCallable,
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


@dataclasses.dataclass
class ResultWorkflows:
    initial_result_workflow: Workflow
    averaging_workflow: Workflow
    forward_output_workflow: Workflow
    base_name: str
    compute_equivalent_before_average: bool = False
    is_single_component_result: bool = False
    mesh_workflow: Optional[Workflow] = None
    principal_workflow: Optional[Workflow] = None
    equivalent_workflow: Optional[Workflow] = None
    norm_workflow: Optional[Workflow] = None
    component_extraction_workflow: Optional[Workflow] = None
    sweeping_phase_workflow: Optional[Workflow] = None


def create_result_workflows(
    base_name: str,
    location: Union[locations, str],
    force_elemental_nodal: bool,
    server,
    mesh: MeshedRegion,
    create_operator_callable: CreateOperatorCallable,
    is_mesh_required: bool,
    has_skin: bool,
    mesh_provider: Any,
    has_principal: bool,
    has_equivalent: bool,
    components_to_extract: list[int],
    should_extract_components: bool,
    has_norm: bool,
    amplitude: bool = False,
    sweeping_phase: Union[float, None] = None,
) -> ResultWorkflows:
    initial_result_wf = create_initial_result_workflow(
        name=base_name,
        location=location,
        force_elemental_nodal=force_elemental_nodal,
        server=server,
        create_operator_callable=create_operator_callable,
        mesh=mesh,
    )

    average_wf = create_averaging_workflow(
        location=location,
        has_skin=has_skin,
        mesh_averaging_needed=force_elemental_nodal,
        create_operator_callable=create_operator_callable,
        server=server,
    )

    forward_output_wf = create_output_workflow(
        create_operator_callable=create_operator_callable, server=server
    )

    result_workflows = ResultWorkflows(
        initial_result_workflow=initial_result_wf,
        averaging_workflow=average_wf,
        forward_output_workflow=forward_output_wf,
        base_name=base_name,
    )

    if is_mesh_required:
        result_workflows.mesh_workflow = create_mesh_workflow(
            mesh_provider=mesh_provider, server=server
        )

    # Add a step to compute principal invariants if result is principal
    if has_principal:
        result_workflows.principal_workflow = create_principal_workflow(
            components_to_extract=components_to_extract,
            create_operator_callable=create_operator_callable,
            server=server,
        )

    # Add a step to compute equivalent if result is equivalent
    if has_equivalent:
        result_workflows.equivalent_workflow = create_equivalent_workflow(
            create_operator_callable=create_operator_callable, server=server
        )
        result_workflows.base_name += "_VM"

        # If a strain result, change the location now
        if has_equivalent and base_name[0] == "E":
            result_workflows.compute_equivalent_before_average = True

    if should_extract_components:
        (
            extract_component_wf,
            base_name,
            result_is_single_component,
        ) = create_extract_component_workflow(
            create_operator_callable=create_operator_callable,
            components_to_extract=components_to_extract,
            base_name=base_name,
            server=server,
        )
        result_workflows.component_extraction_workflow = extract_component_wf
        result_workflows.is_single_component_result = result_is_single_component
        result_workflows.base_name = base_name

    if has_norm:
        norm_wf, base_name = create_norm_workflow(
            create_operator_callable=create_operator_callable,
            base_name=base_name,
            server=server,
        )
        result_workflows.norm_workflow = norm_wf
        result_workflows.is_single_component_result = True
        result_workflows.base_name = base_name

    result_workflows.sweeping_phase_workflow = create_sweeping_phase_workflow(
        create_operator_callable=create_operator_callable,
        server=server,
        amplitude=amplitude,
        sweeping_phase=sweeping_phase,
    )

    return result_workflows


def connect_cyclic_inputs(expand_cyclic, phase_angle_cyclic, result_wf):
    if expand_cyclic is not False:
        # If expand_cyclic is a list
        if isinstance(expand_cyclic, list) and len(expand_cyclic) > 0:
            # If a list of sector numbers, directly connect it to the num_sectors pin
            if all(
                [isinstance(expand_cyclic_i, int) for expand_cyclic_i in expand_cyclic]
            ):
                if any([i < 1 for i in expand_cyclic]):
                    raise ValueError(
                        "Sector selection with 'expand_cyclic' starts at 1."
                    )
                result_wf.connect(
                    _WfNames.cyclic_sectors_to_expand,
                    [i - 1 for i in expand_cyclic],
                )
            # If any is a list, treat it as per stage num_sectors
            elif any(
                [isinstance(expand_cyclic_i, list) for expand_cyclic_i in expand_cyclic]
            ):
                # Create a ScopingsContainer to fill
                sectors_scopings = ScopingsContainer()
                sectors_scopings.labels = ["stage"]
                # For each potential num_sectors, check either an int or a list of ints
                for i, num_sectors_stage_i in enumerate(expand_cyclic):
                    # Prepare num_sectors data
                    if isinstance(num_sectors_stage_i, int):
                        num_sectors_stage_i = [num_sectors_stage_i]
                    elif isinstance(num_sectors_stage_i, list):
                        if not all([isinstance(n, int) for n in num_sectors_stage_i]):
                            raise ValueError(
                                "'expand_cyclic' only accepts lists of int values >= 1."
                            )
                    # num_sectors_stage_i is now a list of int,
                    # add an equivalent Scoping with the correct 'stage' label value
                    if any([i < 1 for i in num_sectors_stage_i]):
                        raise ValueError(
                            "Sector selection with 'expand_cyclic' starts at 1."
                        )
                    sectors_scopings.add_scoping(
                        {"stage": i},
                        Scoping(ids=[i - 1 for i in num_sectors_stage_i]),
                    )
                result_wf.connect(
                    _WfNames.cyclic_sectors_to_expand, inpt=sectors_scopings
                )
        elif not isinstance(expand_cyclic, bool):
            raise ValueError(
                "'expand_cyclic' argument can only be a boolean or a list."
            )
        result_wf.connect(_WfNames.read_cyclic, 3)  # Connect the read_cyclic pin
    else:
        result_wf.connect(_WfNames.read_cyclic, 1)  # Connect the read_cyclic pin
    if phase_angle_cyclic is not None:
        if isinstance(phase_angle_cyclic, int):
            phase_angle_cyclic = float(phase_angle_cyclic)
        if not isinstance(phase_angle_cyclic, float):
            raise ValueError(
                "'phase_angle_cyclic' argument only accepts a single float value."
            )
        result_wf.connect(_WfNames.cyclic_phase, phase_angle_cyclic)
    return result_wf


def connect_initial_results_inputs(
    initial_result_workflow: Workflow,
    selection: Selection,
    data_sources: Any,
    streams_provider: Any,
    expand_cyclic: bool,
    phase_angle_cyclic: Any,
):
    initial_result_workflow.connect_with(
        selection.spatial_selection._selection,
        output_input_names={"scoping": "mesh_scoping"},
    )

    initial_result_workflow.connect_with(
        selection.time_freq_selection._selection,
        output_input_names=("scoping", "time_scoping"),
    )

    # Connect data_sources and streams_container inputs of selection if necessary
    # Todo: was this working before?
    if "streams" in initial_result_workflow.input_names:
        initial_result_workflow.connect("streams", streams_provider)
    if "data_sources" in initial_result_workflow.input_names:
        initial_result_workflow.connect("data_sources", data_sources)

    connect_cyclic_inputs(
        expand_cyclic=expand_cyclic,
        phase_angle_cyclic=phase_angle_cyclic,
        result_wf=initial_result_workflow,
    )


def connect_averaging_eqv_and_principal_workflows(
    result_workflows: ResultWorkflows,
):
    averaging_wf_connections = {
        _WfNames.output_data: _WfNames.input_data,
        _WfNames.skin: _WfNames.skin,
        _WfNames.skin_input_mesh: _WfNames.skin_input_mesh,
    }
    assert not (
        result_workflows.equivalent_workflow is not None
        and result_workflows.principal_workflow is not None
    )

    principal_or_eqv_wf = (
        result_workflows.equivalent_workflow or result_workflows.principal_workflow
    )

    if not result_workflows.compute_equivalent_before_average:
        result_workflows.averaging_workflow.connect_with(
            result_workflows.initial_result_workflow,
            output_input_names=averaging_wf_connections,
        )
        if principal_or_eqv_wf is not None:
            principal_or_eqv_wf.connect_with(
                result_workflows.averaging_workflow,
                output_input_names={_WfNames.output_data: _WfNames.input_data},
            )
            output_wf = principal_or_eqv_wf
        else:
            output_wf = result_workflows.averaging_workflow

    else:
        assert principal_or_eqv_wf is not None
        principal_or_eqv_wf.connect_with(
            result_workflows.initial_result_workflow,
            output_input_names={_WfNames.output_data: _WfNames.input_data},
        )
        result_workflows.averaging_workflow.connect_with(
            principal_or_eqv_wf,
            output_input_names=averaging_wf_connections,
        )
        output_wf = result_workflows.averaging_workflow

    return output_wf


def append_workflow(new_wf: Optional[Workflow], last_wf: Workflow):
    if new_wf is not None:
        new_wf.connect_with(
            last_wf,
            output_input_names={_WfNames.output_data: _WfNames.input_data},
        )
        return new_wf
    else:
        return last_wf
