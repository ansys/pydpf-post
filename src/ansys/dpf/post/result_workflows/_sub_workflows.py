# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from typing import Optional, Union

from ansys.dpf.core import (
    MeshedRegion,
    StreamsContainer,
    Workflow,
    operators,
    shell_layers,
)
from ansys.dpf.gate.common import locations

from ansys.dpf.post.misc import _connect_any
from ansys.dpf.post.result_workflows._utils import _CreateOperatorCallable, _Rescoping
from ansys.dpf.post.selection import SpatialSelection, _WfNames


def _create_averaging_workflow(
    has_skin: bool,
    location: str,
    force_elemental_nodal: bool,
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

    map_to_skin = True
    if not server.meet_version("8.0"):
        # Before 8.0, the skin mapping was not supported
        # for elemental and nodal results (only elemental nodal).
        # In the nodal case the mapping is not needed, but we still
        # call the operator to be consistent.
        # For elemental results the mapping was not working before 8.0.
        map_to_skin = force_elemental_nodal

    if has_skin and map_to_skin:
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
            _connect_any(solid_to_skin_operator.inputs.mesh_scoping, skin_mesh_fwd)
            # To keep for retro-compatibility
        else:
            _connect_any(solid_to_skin_operator.inputs.mesh, skin_mesh_fwd)

        if server.meet_version("8.0"):
            # solid mesh_input only supported for server version
            # 8.0 and up
            average_wf.set_input_name(
                _WfNames.skin_input_mesh, solid_to_skin_operator.inputs.solid_mesh
            )

    if (
        location == locations.nodal or location == locations.elemental
    ) and force_elemental_nodal:
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


def _create_extract_component_workflow(
    create_operator_callable: _CreateOperatorCallable,
    components_to_extract: list[int],
    component_names: list[str],
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
        new_base_name = base_name + f"_{component_names[0]}"
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
    name: str,
    server,
    shell_layer: Optional[shell_layers],
    is_nodal: bool,
    create_operator_callable: _CreateOperatorCallable,
):
    initial_result_workflow = Workflow(server=server)

    initial_result_op = create_operator_callable(name=name)

    initial_result_workflow.set_input_name(_WfNames.mesh, initial_result_op, 7)
    initial_result_workflow.set_input_name(_WfNames.location, initial_result_op, 9)

    initial_result_workflow.add_operator(initial_result_op)

    forward_shell_layer_op = operators.utility.forward(server=server)
    initial_result_workflow.add_operator(forward_shell_layer_op)
    initial_result_workflow.set_input_name(_WfNames.shell_layer, forward_shell_layer_op)

    forward_input_mesh_scoping_op = operators.utility.forward(server=server)
    initial_result_workflow.add_operator(forward_input_mesh_scoping_op)

    # The next section is only needed, because the shell_layer selection does not
    # work for elemental and elemental nodal results.
    # If elemental results are requested with a chosen shell layer,
    # the shell layer is not selected and the results are split into solids
    # and shells. Here, we add an additional shell layer selection and merge_shell_solid
    # operator to manually merge the results.
    # Note that we have to skip this step if the location is nodal, because
    # the resulting location is wrong when the shell layer is selected again manually, after
    # it was already selected by the initial result operator.
    if shell_layer is not None and not is_nodal:
        merge_shell_solid_fields = create_operator_callable(
            name="merge::solid_shell_fields"
        )
        initial_result_workflow.add_operator(merge_shell_solid_fields)
        shell_layer_op = operators.utility.change_shell_layers(server=server)
        shell_layer_op.inputs.merge(True)
        initial_result_workflow.add_operator(shell_layer_op)

        initial_result_workflow.set_output_name(
            _WfNames.output_data, merge_shell_solid_fields, 0
        )
        shell_layer_op.inputs.fields_container(
            initial_result_op.outputs.fields_container
        )
        _connect_any(
            shell_layer_op.inputs.e_shell_layer, forward_shell_layer_op.outputs.any
        )

        merge_shell_solid_fields.inputs.fields_container(
            shell_layer_op.outputs.fields_container_as_fields_container
        )

        # End section for elemental results with shell layer selection
    else:
        initial_result_workflow.set_output_name(
            _WfNames.output_data, initial_result_op, 0
        )

    if hasattr(initial_result_op.inputs, "shell_layer"):
        _connect_any(
            initial_result_op.inputs.shell_layer, forward_shell_layer_op.outputs.any
        )

    initial_result_workflow.set_input_name(
        "time_scoping", initial_result_op.inputs.time_scoping
    )

    initial_result_workflow.set_input_name(
        "mesh_scoping", forward_input_mesh_scoping_op
    )
    _connect_any(
        initial_result_op.inputs.mesh_scoping, forward_input_mesh_scoping_op.outputs.any
    )

    initial_result_workflow.set_output_name(
        _WfNames.result_scoping, forward_input_mesh_scoping_op.outputs.any
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
        sweeping_op.connect(3, "deg")
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


def _enrich_mesh_with_property_fields(
    mesh: MeshedRegion,
    property_names: list[str],
    streams_provider: StreamsContainer,
):
    property_operator = operators.metadata.property_field_provider_by_name(
        server=mesh._server
    )
    property_operator.inputs.streams_container(streams_provider)

    for property_name in property_names:
        # Some of the requested properties might already be part of the mesh
        # property fields
        if property_name not in mesh.available_property_fields:
            property_operator.inputs.property_name(property_name)
            property_field = property_operator.eval()

            # Rescope the property field to the element scoping of the mesh
            # to ensure the split by property operator works correctly
            rescope_op = operators.scoping.rescope_property_field(
                mesh_scoping=mesh.elements.scoping,
                fields=property_field,
                server=mesh._server,
            )

            mesh.set_property_field(
                property_name, rescope_op.outputs.fields_as_property_field()
            )


def _create_split_scope_by_body_workflow(server, body_defining_properties: list[str]):
    split_scope_by_body_wf = Workflow(server=server)
    split_scop_op = operators.scoping.split_on_property_type(server=server)
    split_scope_by_body_wf.add_operator(split_scop_op)
    split_scope_by_body_wf.set_input_name(_WfNames.mesh, split_scop_op.inputs.mesh)
    split_scope_by_body_wf.set_input_name(
        _WfNames.scoping_location, split_scop_op.inputs.requested_location
    )
    split_scope_by_body_wf.set_input_name(
        _WfNames.scoping, split_scop_op.inputs.mesh_scoping
    )

    for idx, property_name in enumerate(body_defining_properties):
        split_scop_op.connect(13 + idx, property_name)
    split_scope_by_body_wf.set_output_name(
        _WfNames.scoping, split_scop_op.outputs.mesh_scoping
    )
    split_scope_by_body_wf.set_output_name(
        _WfNames.result_scoping_by_body, split_scop_op.outputs.mesh_scoping
    )
    return split_scope_by_body_wf


def _create_rescoping_workflow(server, rescoping: _Rescoping):
    selection = SpatialSelection(server=server)

    if rescoping.named_selections is not None:
        selection.select_named_selection(rescoping.named_selections)

    if rescoping.node_ids is not None:
        selection.select_nodes(rescoping.node_ids)

    rescoping_wf = Workflow(server=server)

    transpose_scoping_op = operators.scoping.transpose(server=server)
    rescoping_wf.add_operator(transpose_scoping_op)
    transpose_scoping_op.inputs.requested_location(rescoping.requested_location)
    rescoping_wf.set_input_name(
        _WfNames.mesh, transpose_scoping_op.inputs.meshed_region
    )

    rescoping_op = operators.scoping.rescope_fc(server=server)
    rescoping_wf.add_operator(rescoping_op)
    rescoping_op.inputs.mesh_scoping(
        transpose_scoping_op.outputs.mesh_scoping_as_scoping
    )
    rescoping_wf.set_input_name(
        _WfNames.input_data, rescoping_op.inputs.fields_container
    )
    rescoping_wf.set_input_name(
        _WfNames.scoping, transpose_scoping_op.inputs.mesh_scoping
    )
    rescoping_wf.set_output_name(
        _WfNames.output_data, rescoping_op.outputs.fields_container
    )

    rescoping_wf.connect_with(
        selection._selection, output_input_names={_WfNames.scoping: _WfNames.scoping}
    )

    return rescoping_wf
