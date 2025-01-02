from typing import Any, Optional

from ansys.dpf.core import (
    MeshedRegion,
    Scoping,
    ScopingsContainer,
    Workflow,
    shell_layers,
)

from ansys.dpf.post.result_workflows._build_workflow import ResultWorkflows
from ansys.dpf.post.result_workflows._sub_workflows import (
    _enrich_mesh_with_property_fields,
)
from ansys.dpf.post.result_workflows._utils import AveragingConfig
from ansys.dpf.post.selection import Selection, _WfNames


def _connect_cyclic_inputs(expand_cyclic, phase_angle_cyclic, result_wf: Workflow):
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
    elif _WfNames.read_cyclic in result_wf.input_names:
        result_wf.connect(_WfNames.read_cyclic, 1)  # Connect the read_cyclic pin
    if phase_angle_cyclic is not None:
        if isinstance(phase_angle_cyclic, int):
            phase_angle_cyclic = float(phase_angle_cyclic)
        if not isinstance(phase_angle_cyclic, float):
            raise ValueError(
                "'phase_angle_cyclic' argument only accepts a single float value."
            )
        result_wf.connect(_WfNames.cyclic_phase, phase_angle_cyclic)


def _connect_workflow_inputs(
    initial_result_workflow: Workflow,
    split_by_body_workflow: Optional[Workflow],
    rescoping_workflow: Optional[Workflow],
    force_elemental_nodal: bool,
    location: str,
    selection: Selection,
    expand_cyclic: bool,
    phase_angle_cyclic: Any,
    mesh: MeshedRegion,
    streams_provider: Any,
    data_sources: Any,
    averaging_config: AveragingConfig,
    shell_layer: Optional[shell_layers],
):
    """Connects the inputs of the initial result workflow.

    The initial result workflow is the first workflow in the result workflows chain, which
    extracts the raw results from the data sources.
    """
    selection_wf = selection.spatial_selection._selection

    if selection.spatial_selection.requires_mesh:
        selection_wf.connect(_WfNames.initial_mesh, mesh)

    if averaging_config.average_per_body:
        _enrich_mesh_with_property_fields(
            mesh, averaging_config.body_defining_properties, streams_provider
        )

    if split_by_body_workflow is not None:
        split_by_body_workflow.connect(_WfNames.mesh, mesh)
        if force_elemental_nodal:
            split_by_body_workflow.connect(_WfNames.scoping_location, "ElementalNodal")
        else:
            split_by_body_workflow.connect(_WfNames.scoping_location, location)
        split_by_body_workflow.connect_with(
            selection_wf, output_input_names={_WfNames.scoping: _WfNames.scoping}
        )
        selection_wf = split_by_body_workflow

    initial_result_workflow.connect_with(
        selection_wf,
        output_input_names={"scoping": "mesh_scoping"},
    )

    initial_result_workflow.connect_with(
        selection.time_freq_selection._selection,
        output_input_names=("scoping", "time_scoping"),
    )

    # Connect data_sources and streams_container inputs of selection if necessary
    # Note: streams and data_sources inputs are inherited from the selection_workflow
    # connected above
    if (
        "streams" in initial_result_workflow.input_names
        and streams_provider is not None
    ):
        initial_result_workflow.connect("streams", streams_provider)
    if "data_sources" in initial_result_workflow.input_names:
        initial_result_workflow.connect("data_sources", data_sources)

    _connect_cyclic_inputs(
        expand_cyclic=expand_cyclic,
        phase_angle_cyclic=phase_angle_cyclic,
        result_wf=initial_result_workflow,
    )

    if force_elemental_nodal:
        initial_result_workflow.connect(_WfNames.location, "ElementalNodal")
    elif location:
        initial_result_workflow.connect(_WfNames.location, location)

    initial_result_workflow.connect(_WfNames.mesh, mesh)

    if shell_layer is not None:
        if _WfNames.shell_layer not in initial_result_workflow.input_names:
            raise RuntimeError(
                "The shell_layer input is not supported by this workflow."
            )
        initial_result_workflow.connect(_WfNames.shell_layer, shell_layer.value)

    if rescoping_workflow:
        rescoping_workflow.connect(_WfNames.mesh, mesh)
        if _WfNames.data_sources in rescoping_workflow.input_names:
            rescoping_workflow.connect(_WfNames.data_sources, data_sources)


def _connect_averaging_eqv_and_principal_workflows(
    result_workflows: ResultWorkflows,
):
    """Connects the averaging, equivalent, and principal workflows.

    The order of these workflows depends on result_workflows.compute_equivalent_before_average.
    Only one of equivalent_workflow or principal_workflow can be active at the same time.
    """
    averaging_wf_connections = {
        _WfNames.output_data: _WfNames.input_data,
        _WfNames.skin: _WfNames.skin,
        _WfNames.skin_input_mesh: _WfNames.skin_input_mesh,
    }

    if (
        result_workflows.equivalent_workflow is not None
        and result_workflows.principal_workflow is not None
    ):
        raise AssertionError(
            "The equivalent workflow and principal workflow are both not None."
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
        if principal_or_eqv_wf is None:
            raise AssertionError(
                "The equivalent workflow or principal workflow is None."
            )

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
