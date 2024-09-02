from typing import Optional, Protocol

from ansys.dpf.core import Operator, Workflow

from ansys.dpf.post.selection import _WfNames


class _CreateOperatorCallable(Protocol):
    # Callable to create an operator with a given name.
    # This usually corresponds to model.operator
    def __call__(self, name: str) -> Operator:
        ...


def _append_workflows(workflows: list[Workflow], current_output_workflow: Workflow):
    # Append multiple workflows to last_wf. The new
    # workflows must have an _WfNames.input_data and _WfNames.output_data pin. The
    # current_output_workflow must have an _WfNames.output_data pin.
    # Returns the last appended workflow.
    # Workflows are not appended if they are None.
    for workflow in workflows:
        current_output_workflow = _append_workflow(
            new_wf=workflow, last_wf=current_output_workflow
        )
    return current_output_workflow


def _append_workflow(new_wf: Optional[Workflow], last_wf: Workflow):
    # Append a single workflow to last_wf. The new
    # workflow must have an _WfNames.input_data pin and the last_wf
    # must have an _WfNames.output_data pin.
    # Returns the appended workflow if it was not None, otherwise returns last_wf.
    if new_wf is None:
        return last_wf

    assert _WfNames.input_data in new_wf.input_names
    assert _WfNames.output_data in new_wf.output_names
    assert _WfNames.output_data in last_wf.output_names
    new_wf.connect_with(
        last_wf,
        output_input_names={_WfNames.output_data: _WfNames.input_data},
    )
    return new_wf
