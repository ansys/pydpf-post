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

import dataclasses
from typing import Optional, Protocol

from ansys.dpf.core import Operator, Workflow
from ansys.dpf.core.available_result import AvailableResult

from ansys.dpf.post.selection import _WfNames


class _CreateOperatorCallable(Protocol):
    # Callable to create an operator with a given name.
    # This usually corresponds to model.operator
    def __call__(self, name: str) -> Operator:
        ...


class _Rescoping:
    # Defines a rescoping that needs to be performed at the end
    # of the results workflow. This is needed, because
    # the scoping sometimes needs to be broadened when force_elemental_nodal is
    # True.
    def __init__(
        self,
        requested_location: str,
        named_selections: Optional[list[str]] = None,
        node_ids: Optional[list[int]] = None,
    ):
        if named_selections is not None and node_ids is not None:
            raise ValueError(
                "Arguments named_selections and node_ids are mutually exclusive"
            )
        if named_selections is None and node_ids is None:
            raise ValueError(
                "At least one of named_selections and node_ids must be provided"
            )
        self._node_ids = node_ids
        self._named_selections = named_selections
        self._requested_location = requested_location

    @property
    def node_ids(self):
        return self._node_ids

    @property
    def named_selections(self):
        return self._named_selections

    @property
    def requested_location(self):
        return self._requested_location


@dataclasses.dataclass
class AveragingConfig:
    """Configuration for averaging of results."""

    # List of properties that define a body. The mesh is split by these properties to
    # get the bodies.
    body_defining_properties: Optional[list[str]] = None
    # If True, the results are averaged per body. The bodies are determined
    # by the body_defining_properties.
    average_per_body: bool = False


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

    if _WfNames.input_data not in new_wf.input_names:
        raise AssertionError(
            f"Workflow {new_wf} must have an input pin {_WfNames.input_data}"
        )
    if _WfNames.output_data not in new_wf.output_names:
        raise AssertionError(
            f"Workflow {new_wf} must have an output pin {_WfNames.output_data}"
        )
    if _WfNames.output_data not in last_wf.output_names:
        raise AssertionError(
            f"Workflow {last_wf} must have an output pin {_WfNames.output_data}"
        )

    new_wf.connect_with(
        last_wf,
        output_input_names={_WfNames.output_data: _WfNames.input_data},
    )
    return new_wf


def _get_native_location(
    available_results: list[AvailableResult], base_name: str
) -> str | None:
    """Get the native location of a result from its base name."""
    res = next((r for r in available_results if r.operator_name == base_name), None)
    native_location = None
    if res is not None:
        native_location = res.native_location

    return native_location
