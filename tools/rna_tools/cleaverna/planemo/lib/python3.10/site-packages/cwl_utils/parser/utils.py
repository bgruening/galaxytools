"""CWL parser utility functions."""

import copy
import logging
import os
from collections.abc import MutableSequence
from pathlib import Path
from types import ModuleType
from typing import Any, Optional, Union, cast
from urllib.parse import unquote_plus, urlparse

from schema_salad.exceptions import ValidationException
from schema_salad.sourceline import SourceLine, strip_dup_lineno
from schema_salad.utils import json_dumps, yaml_no_ts

import cwl_utils
import cwl_utils.parser

from . import (
    LoadingOptions,
    Process,
    Workflow,
    WorkflowStep,
    WorkflowStepInput,
    cwl_v1_0,
    cwl_v1_0_utils,
    cwl_v1_1,
    cwl_v1_1_utils,
    cwl_v1_2,
    cwl_v1_2_utils,
)

_logger = logging.getLogger("cwl_utils")


def convert_stdstreams_to_files(process: Process) -> None:
    """Convert stdin, stdout and stderr type shortcuts to files."""
    if isinstance(process, cwl_v1_0.CommandLineTool):
        cwl_v1_0_utils.convert_stdstreams_to_files(process)
    elif isinstance(process, cwl_v1_1.CommandLineTool):
        cwl_v1_1_utils.convert_stdstreams_to_files(process)
    elif isinstance(process, cwl_v1_2.CommandLineTool):
        cwl_v1_2_utils.convert_stdstreams_to_files(process)


def load_inputfile_by_uri(
    version: str,
    path: Union[str, Path],
    loadingOptions: Optional[LoadingOptions] = None,
) -> Any:
    """Load a CWL input file from a URI or a path."""
    if isinstance(path, str):
        uri = urlparse(path)
        if not uri.scheme or uri.scheme == "file":
            real_path = Path(unquote_plus(uri.path)).resolve().as_uri()
        else:
            real_path = path
    else:
        real_path = path.resolve().as_uri()

    if version is None:
        raise ValidationException("could not get the cwlVersion")

    baseuri = str(real_path)

    if loadingOptions is None:
        if version == "v1.0":
            loadingOptions = cwl_v1_0.LoadingOptions(fileuri=baseuri)
        elif version == "v1.1":
            loadingOptions = cwl_v1_1.LoadingOptions(fileuri=baseuri)
        elif version == "v1.2":
            loadingOptions = cwl_v1_2.LoadingOptions(fileuri=baseuri)
        else:
            raise ValidationException(
                f"Version error. Did not recognise {version} as a CWL version"
            )

    doc = loadingOptions.fetcher.fetch_text(real_path)
    return load_inputfile_by_string(version, doc, baseuri, loadingOptions)


def load_inputfile(
    version: str,
    doc: Any,
    baseuri: Optional[str] = None,
    loadingOptions: Optional[LoadingOptions] = None,
) -> Any:
    """Load a CWL input file from a serialized YAML string or a YAML object."""
    if baseuri is None:
        baseuri = cwl_v1_0.file_uri(os.getcwd()) + "/"
    if isinstance(doc, str):
        return load_inputfile_by_string(version, doc, baseuri, loadingOptions)
    return load_inputfile_by_yaml(version, doc, baseuri, loadingOptions)


def load_inputfile_by_string(
    version: str,
    string: str,
    uri: str,
    loadingOptions: Optional[LoadingOptions] = None,
) -> Any:
    """Load a CWL input file from a serialized YAML string."""
    yaml = yaml_no_ts()
    result = yaml.load(string)
    return load_inputfile_by_yaml(version, result, uri, loadingOptions)


def load_inputfile_by_yaml(
    version: str,
    yaml: Any,
    uri: str,
    loadingOptions: Optional[LoadingOptions] = None,
) -> Any:
    """Load a CWL input file from a YAML object."""
    if version == "v1.0":
        result = cwl_v1_0_utils.load_inputfile_by_yaml(
            yaml, uri, cast(Optional[cwl_v1_0.LoadingOptions], loadingOptions)
        )
    elif version == "v1.1":
        result = cwl_v1_1_utils.load_inputfile_by_yaml(
            yaml, uri, cast(Optional[cwl_v1_1.LoadingOptions], loadingOptions)
        )
    elif version == "v1.2":
        result = cwl_v1_2_utils.load_inputfile_by_yaml(
            yaml, uri, cast(Optional[cwl_v1_2.LoadingOptions], loadingOptions)
        )
    elif version is None:
        raise ValidationException("could not get the cwlVersion")
    else:
        raise ValidationException(
            f"Version error. Did not recognise {version} as a CWL version"
        )

    return result


def load_step(
    step: cwl_utils.parser.WorkflowStep,
) -> Process:
    if isinstance(step.run, str):
        step_run = cwl_utils.parser.load_document_by_uri(
            path=step.loadingOptions.fetcher.urljoin(
                base_url=cast(str, step.loadingOptions.fileuri),
                url=step.run,
            ),
            loadingOptions=step.loadingOptions,
        )
        return cast(Process, step_run)
    else:
        return cast(Process, copy.deepcopy(step.run))


def static_checker(workflow: cwl_utils.parser.Workflow) -> None:
    """Check if all source and sink types of a workflow are compatible before run time."""
    step_inputs = []
    step_outputs = []
    type_dict = {}
    param_to_step = {}
    for step in workflow.steps:
        if step.in_ is not None:
            step_inputs.extend(step.in_)
            param_to_step.update({s.id: step for s in step.in_})
            type_dict.update(
                {
                    cast(str, s.id): type_for_step_input(
                        step, s, cast(str, workflow.cwlVersion)
                    )
                    for s in step.in_
                }
            )
        if step.out is not None:
            # FIXME: the correct behaviour here would be to create WorkflowStepOutput directly at load time
            if workflow.cwlVersion == "v1.0":
                step_outs = [
                    cwl_v1_0.WorkflowStepOutput(s) if isinstance(s, str) else s
                    for s in step.out
                ]
            elif workflow.cwlVersion == "v1.1":
                step_outs = [
                    cwl_v1_1.WorkflowStepOutput(s) if isinstance(s, str) else s
                    for s in step.out
                ]
            elif workflow.cwlVersion == "v1.2":
                step_outs = [
                    cwl_v1_2.WorkflowStepOutput(s) if isinstance(s, str) else s
                    for s in step.out
                ]
            else:
                raise Exception(f"Unsupported CWL version {workflow.cwlVersion}")
            step_outputs.extend(step_outs)
            param_to_step.update({s.id: step for s in step_outs})
            type_dict.update(
                {
                    s.id: type_for_step_output(step, s.id, workflow.cwlVersion)
                    for s in step_outs
                }
            )
    src_dict = {
        **{param.id: param for param in workflow.inputs},
        **{param.id: param for param in step_outputs},
    }
    type_dict = {
        **type_dict,
        **{param.id: param.type_ for param in workflow.inputs},
        **{param.id: param.type_ for param in workflow.outputs},
    }

    parser: ModuleType
    step_inputs_val: dict[str, Any]
    workflow_outputs_val: dict[str, Any]
    if workflow.cwlVersion == "v1.0":
        parser = cwl_v1_0
        step_inputs_val = cwl_v1_0_utils.check_all_types(
            src_dict, step_inputs, type_dict
        )
        workflow_outputs_val = cwl_v1_0_utils.check_all_types(
            src_dict, workflow.outputs, type_dict
        )
    elif workflow.cwlVersion == "v1.1":
        parser = cwl_v1_1
        step_inputs_val = cwl_v1_1_utils.check_all_types(
            src_dict, step_inputs, type_dict
        )
        workflow_outputs_val = cwl_v1_1_utils.check_all_types(
            src_dict, workflow.outputs, type_dict
        )
    elif workflow.cwlVersion == "v1.2":
        parser = cwl_v1_2
        step_inputs_val = cwl_v1_2_utils.check_all_types(
            src_dict, step_inputs, param_to_step, type_dict
        )
        workflow_outputs_val = cwl_v1_2_utils.check_all_types(
            src_dict, workflow.outputs, param_to_step, type_dict
        )
    else:
        raise Exception(f"Unsupported CWL version {workflow.cwlVersion}")

    warnings = step_inputs_val["warning"] + workflow_outputs_val["warning"]
    exceptions = step_inputs_val["exception"] + workflow_outputs_val["exception"]

    warning_msgs = []
    exception_msgs = []
    for warning in warnings:
        src = warning.src
        sink = warning.sink
        linkMerge = warning.linkMerge
        msg = (
            SourceLine(src, "type").makeError(
                "Source '%s' of type %s may be incompatible"
                % (
                    parser.shortname(src.id),
                    json_dumps(parser.save(type_dict[src.id])),
                )
            )
            + "\n"
            + SourceLine(sink, "type").makeError(
                "  with sink '%s' of type %s"
                % (
                    parser.shortname(sink.id),
                    json_dumps(parser.save(type_dict[sink.id])),
                )
            )
        )
        if linkMerge is not None:
            msg += "\n" + SourceLine(sink).makeError(
                "  source has linkMerge method %s" % linkMerge
            )

        if warning.message is not None:
            msg += "\n" + SourceLine(sink).makeError("  " + warning.message)

        if msg:
            warning_msgs.append(msg)

    for exception in exceptions:
        src = exception.src
        sink = exception.sink
        linkMerge = exception.linkMerge
        extra_message = exception.message

        msg = (
            SourceLine(src, "type").makeError(
                "Source '%s' of type %s is incompatible"
                % (parser.shortname(src.id), json_dumps(parser.save(type_dict[src.id])))
            )
            + "\n"
            + SourceLine(sink, "type").makeError(
                "  with sink '%s' of type %s"
                % (
                    parser.shortname(sink.id),
                    json_dumps(parser.save(type_dict[sink.id])),
                )
            )
        )
        if extra_message is not None:
            msg += "\n" + SourceLine(sink).makeError("  " + extra_message)

        if linkMerge is not None:
            msg += "\n" + SourceLine(sink).makeError(
                "  source has linkMerge method %s" % linkMerge
            )
        exception_msgs.append(msg)

    for sink in step_inputs:
        if (
            "null" != type_dict[sink.id]
            and not (
                isinstance(type_dict[sink.id], MutableSequence)
                and "null" in type_dict[sink.id]
            )
            and getattr(sink, "source", None) is None
            and getattr(sink, "default", None) is None
            and getattr(sink, "valueFrom", None) is None
        ):
            msg = SourceLine(sink).makeError(
                "Required parameter '%s' does not have source, default, or valueFrom expression"
                % parser.shortname(sink.id)
            )
            exception_msgs.append(msg)

    all_warning_msg = strip_dup_lineno("\n".join(warning_msgs))
    all_exception_msg = strip_dup_lineno("\n" + "\n".join(exception_msgs))

    if all_warning_msg:
        _logger.warning("Workflow checker warning:\n%s", all_warning_msg)
    if exceptions:
        raise ValidationException(all_exception_msg)


def type_for_source(
    process: Process,
    sourcenames: Union[str, list[str]],
    parent: Optional[Workflow] = None,
    linkMerge: Optional[str] = None,
    pickValue: Optional[str] = None,
) -> Any:
    """Determine the type for the given sourcenames."""
    if process.cwlVersion == "v1.0":
        return cwl_v1_0_utils.type_for_source(
            cast(
                Union[
                    cwl_v1_0.CommandLineTool,
                    cwl_v1_0.Workflow,
                    cwl_v1_0.ExpressionTool,
                ],
                process,
            ),
            sourcenames,
            cast(Optional[cwl_v1_0.Workflow], parent),
            linkMerge,
        )
    elif process.cwlVersion == "v1.1":
        return cwl_v1_1_utils.type_for_source(
            cast(
                Union[
                    cwl_v1_1.CommandLineTool,
                    cwl_v1_1.Workflow,
                    cwl_v1_1.ExpressionTool,
                ],
                process,
            ),
            sourcenames,
            cast(Optional[cwl_v1_1.Workflow], parent),
            linkMerge,
        )
    elif process.cwlVersion == "v1.2":
        return cwl_v1_2_utils.type_for_source(
            cast(
                Union[
                    cwl_v1_2.CommandLineTool,
                    cwl_v1_2.Workflow,
                    cwl_v1_2.ExpressionTool,
                ],
                process,
            ),
            sourcenames,
            cast(Optional[cwl_v1_2.Workflow], parent),
            linkMerge,
            pickValue,
        )
    elif process.cwlVersion is None:
        raise ValidationException("could not get the cwlVersion")
    else:
        raise ValidationException(
            f"Version error. Did not recognise {process.cwlVersion} as a CWL version"
        )


def type_for_step_input(
    step: WorkflowStep, in_: WorkflowStepInput, cwlVersion: str
) -> Any:
    """Determine the type for the given step output."""
    if cwlVersion == "v1.0":
        return cwl_v1_0_utils.type_for_step_input(
            cast(cwl_v1_0.WorkflowStep, step), cast(cwl_v1_0.WorkflowStepInput, in_)
        )
    elif cwlVersion == "v1.1":
        return cwl_v1_1_utils.type_for_step_input(
            cast(cwl_v1_1.WorkflowStep, step), cast(cwl_v1_1.WorkflowStepInput, in_)
        )
    elif cwlVersion == "v1.2":
        return cwl_v1_2_utils.type_for_step_input(
            cast(cwl_v1_2.WorkflowStep, step), cast(cwl_v1_2.WorkflowStepInput, in_)
        )


def type_for_step_output(step: WorkflowStep, sourcename: str, cwlVersion: str) -> Any:
    """Determine the type for the given step output."""
    if cwlVersion == "v1.0":
        return cwl_v1_0_utils.type_for_step_output(
            cast(cwl_v1_0.WorkflowStep, step), sourcename
        )
    elif cwlVersion == "v1.1":
        return cwl_v1_1_utils.type_for_step_output(
            cast(cwl_v1_1.WorkflowStep, step), sourcename
        )
    elif cwlVersion == "v1.2":
        return cwl_v1_2_utils.type_for_step_output(
            cast(cwl_v1_2.WorkflowStep, step), sourcename
        )


def param_for_source_id(
    process: Union[
        cwl_utils.parser.CommandLineTool,
        cwl_utils.parser.Workflow,
        cwl_utils.parser.ExpressionTool,
    ],
    sourcenames: Union[str, list[str]],
    parent: Optional[cwl_utils.parser.Workflow] = None,
    scatter_context: Optional[list[Optional[tuple[int, str]]]] = None,
) -> Union[
    Union[
        list[cwl_utils.parser.cwl_v1_0.InputParameter],
        cwl_utils.parser.cwl_v1_0.InputParameter,
    ],
    Union[
        list[cwl_utils.parser.cwl_v1_1.WorkflowInputParameter],
        cwl_utils.parser.cwl_v1_1.WorkflowInputParameter,
    ],
    Union[
        list[cwl_utils.parser.cwl_v1_2.WorkflowInputParameter],
        cwl_utils.parser.cwl_v1_2.WorkflowInputParameter,
    ],
]:
    if process.cwlVersion == "v1.0":
        return cwl_utils.parser.cwl_v1_0_utils.param_for_source_id(
            cast(
                Union[
                    cwl_utils.parser.cwl_v1_0.CommandLineTool,
                    cwl_utils.parser.cwl_v1_0.Workflow,
                    cwl_utils.parser.cwl_v1_0.ExpressionTool,
                ],
                process,
            ),
            sourcenames,
            cast(cwl_utils.parser.cwl_v1_0.Workflow, parent),
            scatter_context,
        )
    elif process.cwlVersion == "v1.1":
        return cwl_utils.parser.cwl_v1_1_utils.param_for_source_id(
            cast(
                Union[
                    cwl_utils.parser.cwl_v1_1.CommandLineTool,
                    cwl_utils.parser.cwl_v1_1.Workflow,
                    cwl_utils.parser.cwl_v1_1.ExpressionTool,
                ],
                process,
            ),
            sourcenames,
            cast(cwl_utils.parser.cwl_v1_1.Workflow, parent),
            scatter_context,
        )
    elif process.cwlVersion == "v1.2":
        return cwl_utils.parser.cwl_v1_2_utils.param_for_source_id(
            cast(
                Union[
                    cwl_utils.parser.cwl_v1_2.CommandLineTool,
                    cwl_utils.parser.cwl_v1_2.Workflow,
                    cwl_utils.parser.cwl_v1_2.ExpressionTool,
                ],
                process,
            ),
            sourcenames,
            cast(cwl_utils.parser.cwl_v1_2.Workflow, parent),
            scatter_context,
        )
    elif process.cwlVersion is None:
        raise ValidationException("could not get the cwlVersion")
    else:
        raise ValidationException(
            f"Version error. Did not recognise {process.cwlVersion} as a CWL version"
        )
