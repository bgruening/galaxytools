# SPDX-License-Identifier: Apache-2.0
import hashlib
import logging
import os
from collections import namedtuple
from collections.abc import MutableMapping, MutableSequence
from io import StringIO
from typing import IO, Any, Optional, Union, cast
from urllib.parse import urldefrag

from schema_salad.exceptions import ValidationException
from schema_salad.sourceline import SourceLine, add_lc_filename
from schema_salad.utils import aslist, json_dumps, yaml_no_ts

import cwl_utils.parser
import cwl_utils.parser.cwl_v1_0 as cwl
import cwl_utils.parser.utils
from cwl_utils.errors import WorkflowException
from cwl_utils.utils import yaml_dumps

CONTENT_LIMIT: int = 64 * 1024

_logger = logging.getLogger("cwl_utils")

SrcSink = namedtuple("SrcSink", ["src", "sink", "linkMerge", "message"])


def _compare_records(
    src: cwl.RecordSchema, sink: cwl.RecordSchema, strict: bool = False
) -> bool:
    """
    Compare two records, ensuring they have compatible fields.

    This handles normalizing record names, which will be relative to workflow
    step, so that they can be compared.
    """
    srcfields = {cwl.shortname(field.name): field.type_ for field in (src.fields or {})}
    sinkfields = {
        cwl.shortname(field.name): field.type_ for field in (sink.fields or {})
    }
    for key in sinkfields.keys():
        if (
            not can_assign_src_to_sink(
                srcfields.get(key, "null"), sinkfields.get(key, "null"), strict
            )
            and sinkfields.get(key) is not None
        ):
            _logger.info(
                "Record comparison failure for %s and %s\n"
                "Did not match fields for %s: %s and %s",
                cast(
                    Union[cwl.InputRecordSchema, cwl.CommandOutputRecordSchema], src
                ).name,
                cast(
                    Union[cwl.InputRecordSchema, cwl.CommandOutputRecordSchema], sink
                ).name,
                key,
                srcfields.get(key),
                sinkfields.get(key),
            )
            return False
    return True


def _compare_type(type1: Any, type2: Any) -> bool:
    if isinstance(type1, cwl.ArraySchema) and isinstance(type2, cwl.ArraySchema):
        return _compare_type(type1.items, type2.items)
    elif isinstance(type1, cwl.RecordSchema) and isinstance(type2, cwl.RecordSchema):
        fields1 = {
            cwl.shortname(field.name): field.type_ for field in (type1.fields or {})
        }
        fields2 = {
            cwl.shortname(field.name): field.type_ for field in (type2.fields or {})
        }
        if fields1.keys() != fields2.keys():
            return False
        return all(_compare_type(fields1[k], fields2[k]) for k in fields1.keys())
    elif isinstance(type1, MutableSequence) and isinstance(type2, MutableSequence):
        if len(type1) != len(type2):
            return False
        for t1 in type1:
            if not any(_compare_type(t1, t2) for t2 in type2):
                return False
        return True
    else:
        return bool(type1 == type2)


def _inputfile_load(
    doc: Union[str, MutableMapping[str, Any], MutableSequence[Any]],
    baseuri: str,
    loadingOptions: cwl.LoadingOptions,
    addl_metadata_fields: Optional[MutableSequence[str]] = None,
) -> tuple[Any, cwl.LoadingOptions]:
    loader = cwl.CWLInputFileLoader
    if isinstance(doc, str):
        url = loadingOptions.fetcher.urljoin(baseuri, doc)
        if url in loadingOptions.idx:
            return loadingOptions.idx[url]
        doc_url, frg = urldefrag(url)
        text = loadingOptions.fetcher.fetch_text(doc_url)
        textIO = StringIO(text)
        textIO.name = str(doc_url)
        yaml = yaml_no_ts()
        result = yaml.load(textIO)
        add_lc_filename(result, doc_url)
        loadingOptions = cwl.LoadingOptions(copyfrom=loadingOptions, fileuri=doc_url)
        _inputfile_load(
            result,
            doc_url,
            loadingOptions,
        )
        return loadingOptions.idx[url]

    if isinstance(doc, MutableMapping):
        addl_metadata = {}
        if addl_metadata_fields is not None:
            for mf in addl_metadata_fields:
                if mf in doc:
                    addl_metadata[mf] = doc[mf]

        loadingOptions = cwl.LoadingOptions(
            copyfrom=loadingOptions,
            baseuri=baseuri,
            addl_metadata=addl_metadata,
        )

        loadingOptions.idx[baseuri] = (
            loader.load(doc, baseuri, loadingOptions, docRoot=baseuri),
            loadingOptions,
        )

        return loadingOptions.idx[baseuri]

    if isinstance(doc, MutableSequence):
        loadingOptions.idx[baseuri] = (
            loader.load(doc, baseuri, loadingOptions),
            loadingOptions,
        )
        return loadingOptions.idx[baseuri]

    raise ValidationException(
        "Expected URI string, MutableMapping or MutableSequence, got %s" % type(doc)
    )


def can_assign_src_to_sink(src: Any, sink: Any, strict: bool = False) -> bool:
    """
    Check for identical type specifications, ignoring extra keys like inputBinding.

    src: admissible source types
    sink: admissible sink types

    In non-strict comparison, at least one source type must match one sink type,
       except for 'null'.
    In strict comparison, all source types must match at least one sink type.
    """
    if src == "Any" or sink == "Any":
        return True
    if isinstance(src, cwl.ArraySchema) and isinstance(sink, cwl.ArraySchema):
        return can_assign_src_to_sink(src.items, sink.items, strict)
    if isinstance(src, cwl.RecordSchema) and isinstance(sink, cwl.RecordSchema):
        return _compare_records(src, sink, strict)
    if isinstance(src, MutableSequence):
        if strict:
            for this_src in src:
                if not can_assign_src_to_sink(this_src, sink):
                    return False
            return True
        for this_src in src:
            if this_src != "null" and can_assign_src_to_sink(this_src, sink):
                return True
        return False
    if isinstance(sink, MutableSequence):
        for this_sink in sink:
            if can_assign_src_to_sink(src, this_sink):
                return True
        return False
    return bool(src == sink)


def check_all_types(
    src_dict: dict[str, Any],
    sinks: MutableSequence[Union[cwl.WorkflowStepInput, cwl.WorkflowOutputParameter]],
    type_dict: dict[str, Any],
) -> dict[str, list[SrcSink]]:
    """Given a list of sinks, check if their types match with the types of their sources."""
    validation: dict[str, list[SrcSink]] = {"warning": [], "exception": []}
    for sink in sinks:
        if isinstance(sink, cwl.WorkflowOutputParameter):
            sourceName = "outputSource"
            sourceField = sink.outputSource
        elif isinstance(sink, cwl.WorkflowStepInput):
            sourceName = "source"
            sourceField = sink.source
        else:
            continue
        if sourceField is not None:
            if isinstance(sourceField, MutableSequence):
                linkMerge = sink.linkMerge or (
                    "merge_nested" if len(sourceField) > 1 else None
                )
                srcs_of_sink = []
                for parm_id in sourceField:
                    srcs_of_sink += [src_dict[parm_id]]
            else:
                parm_id = cast(str, sourceField)
                if parm_id not in src_dict:
                    raise SourceLine(sink, sourceName, ValidationException).makeError(
                        f"{sourceName} not found: {parm_id}"
                    )
                srcs_of_sink = [src_dict[parm_id]]
                linkMerge = None
            for src in srcs_of_sink:
                check_result = check_types(
                    type_dict[cast(str, src.id)],
                    type_dict[sink.id],
                    linkMerge,
                    getattr(sink, "valueFrom", None),
                )
                if check_result == "warning":
                    validation["warning"].append(SrcSink(src, sink, linkMerge, None))
                elif check_result == "exception":
                    validation["exception"].append(SrcSink(src, sink, linkMerge, None))
    return validation


def check_types(
    srctype: Any,
    sinktype: Any,
    linkMerge: Optional[str],
    valueFrom: Optional[str] = None,
) -> str:
    """
    Check if the source and sink types are correct.

    Acceptable types are "pass", "warning", or "exception".
    """
    if valueFrom is not None:
        return "pass"
    if linkMerge is None:
        if can_assign_src_to_sink(srctype, sinktype, strict=True):
            return "pass"
        if can_assign_src_to_sink(srctype, sinktype, strict=False):
            return "warning"
        return "exception"
    if linkMerge == "merge_nested":
        return check_types(
            cwl.ArraySchema(items=srctype, type_="array"), sinktype, None, None
        )
    if linkMerge == "merge_flattened":
        return check_types(merge_flatten_type(srctype), sinktype, None, None)
    raise ValidationException(f"Invalid value {linkMerge} for linkMerge field.")


def content_limit_respected_read_bytes(f: IO[bytes]) -> bytes:
    """
    Read file content up to 64 kB as a byte array.

    Truncate content for larger files.
    """
    return f.read(CONTENT_LIMIT)


def content_limit_respected_read(f: IO[bytes]) -> str:
    """
    Read file content up to 64 kB as an utf-8 encoded string.

    Truncate content for larger files.
    """
    return content_limit_respected_read_bytes(f).decode("utf-8")


def convert_stdstreams_to_files(clt: cwl.CommandLineTool) -> None:
    """Convert stdout and stderr type shortcuts to files."""
    for out in clt.outputs:
        if out.type_ == "stdout":
            if out.outputBinding is not None:
                raise ValidationException(
                    "Not allowed to specify outputBinding when using stdout shortcut."
                )
            if clt.stdout is None:
                clt.stdout = str(
                    hashlib.sha1(  # nosec
                        json_dumps(clt.save(), sort_keys=True).encode("utf-8")
                    ).hexdigest()
                )
            out.type_ = "File"
            out.outputBinding = cwl.CommandOutputBinding(glob=clt.stdout)
        elif out.type_ == "stderr":
            if out.outputBinding is not None:
                raise ValidationException(
                    "Not allowed to specify outputBinding when using stderr shortcut."
                )
            if clt.stderr is None:
                clt.stderr = str(
                    hashlib.sha1(  # nosec
                        json_dumps(clt.save(), sort_keys=True).encode("utf-8")
                    ).hexdigest()
                )
            out.type_ = "File"
            out.outputBinding = cwl.CommandOutputBinding(glob=clt.stderr)


def load_inputfile(
    doc: Any,
    baseuri: Optional[str] = None,
    loadingOptions: Optional[cwl.LoadingOptions] = None,
) -> Any:
    """Load a CWL v1.0 input file from a serialized YAML string or a YAML object."""
    if baseuri is None:
        baseuri = cwl.file_uri(os.getcwd()) + "/"
    if loadingOptions is None:
        loadingOptions = cwl.LoadingOptions()

    result, metadata = _inputfile_load(
        doc,
        baseuri,
        loadingOptions,
    )
    return result


def load_inputfile_by_string(
    string: Any,
    uri: str,
    loadingOptions: Optional[cwl.LoadingOptions] = None,
) -> Any:
    """Load a CWL v1.0 input file from a serialized YAML string."""
    yaml = yaml_no_ts()
    result = yaml.load(string)
    add_lc_filename(result, uri)

    if loadingOptions is None:
        loadingOptions = cwl.LoadingOptions(fileuri=uri)

    result, metadata = _inputfile_load(
        result,
        uri,
        loadingOptions,
    )
    return result


def load_inputfile_by_yaml(
    yaml: Any,
    uri: str,
    loadingOptions: Optional[cwl.LoadingOptions] = None,
) -> Any:
    """Load a CWL v1.0 input file from a YAML object."""
    add_lc_filename(yaml, uri)

    if loadingOptions is None:
        loadingOptions = cwl.LoadingOptions(fileuri=uri)

    result, metadata = _inputfile_load(
        yaml,
        uri,
        loadingOptions,
    )
    return result


def merge_flatten_type(src: Any) -> Any:
    """Return the merge flattened type of the source type."""
    if isinstance(src, MutableSequence):
        return [merge_flatten_type(t) for t in src]
    if isinstance(src, cwl.ArraySchema):
        return src
    return cwl.ArraySchema(type_="array", items=src)


def type_for_step_input(
    step: cwl.WorkflowStep,
    in_: cwl.WorkflowStepInput,
) -> Any:
    """Determine the type for the given step input."""
    if in_.valueFrom is not None:
        return "Any"
    step_run = cwl_utils.parser.utils.load_step(step)
    cwl_utils.parser.utils.convert_stdstreams_to_files(step_run)
    if step_run and step_run.inputs:
        for step_input in step_run.inputs:
            if cast(str, step_input.id).split("#")[-1] == in_.id.split("#")[-1]:
                input_type = step_input.type_
                if step.scatter is not None and in_.id in aslist(step.scatter):
                    input_type = cwl.ArraySchema(items=input_type, type_="array")
                return input_type
    return "Any"


def type_for_step_output(
    step: cwl.WorkflowStep,
    sourcename: str,
) -> Any:
    """Determine the type for the given step output."""
    step_run = cwl_utils.parser.utils.load_step(step)
    cwl_utils.parser.utils.convert_stdstreams_to_files(step_run)
    if step_run and step_run.outputs:
        for step_output in step_run.outputs:
            if (
                step_output.id.split("#")[-1].split("/")[-1]
                == sourcename.split("#")[-1].split("/")[-1]
            ):
                output_type = step_output.type_
                if step.scatter is not None:
                    if step.scatterMethod == "nested_crossproduct":
                        for _ in range(len(aslist(step.scatter))):
                            output_type = cwl.ArraySchema(
                                items=output_type, type_="array"
                            )
                    else:
                        output_type = cwl.ArraySchema(items=output_type, type_="array")
                return output_type
    raise ValidationException(
        "param {} not found in {}.".format(
            sourcename,
            yaml_dumps(cwl.save(step)),
        )
    )


def type_for_source(
    process: Union[cwl.CommandLineTool, cwl.Workflow, cwl.ExpressionTool],
    sourcenames: Union[str, list[str]],
    parent: Optional[cwl.Workflow] = None,
    linkMerge: Optional[str] = None,
) -> Any:
    """Determine the type for the given sourcenames."""
    scatter_context: list[Optional[tuple[int, str]]] = []
    params = param_for_source_id(process, sourcenames, parent, scatter_context)
    if not isinstance(params, list):
        new_type = params.type_
        if scatter_context[0] is not None:
            if scatter_context[0][1] == "nested_crossproduct":
                for _ in range(scatter_context[0][0]):
                    new_type = cwl.ArraySchema(items=new_type, type_="array")
            else:
                new_type = cwl.ArraySchema(items=new_type, type_="array")
        if linkMerge == "merge_nested":
            new_type = cwl.ArraySchema(items=new_type, type_="array")
        elif linkMerge == "merge_flattened":
            new_type = merge_flatten_type(new_type)
        return new_type
    new_type = []
    for p, sc in zip(params, scatter_context):
        if isinstance(p, str) and not any(_compare_type(t, p) for t in new_type):
            cur_type = p
        elif hasattr(p, "type_") and not any(
            _compare_type(t, p.type_) for t in new_type
        ):
            cur_type = p.type_
        else:
            cur_type = None
        if cur_type is not None:
            if sc is not None:
                if sc[1] == "nested_crossproduct":
                    for _ in range(sc[0]):
                        cur_type = cwl.ArraySchema(items=cur_type, type_="array")
                else:
                    cur_type = cwl.ArraySchema(items=cur_type, type_="array")
            new_type.append(cur_type)
    if len(new_type) == 1:
        new_type = new_type[0]
    if linkMerge == "merge_nested":
        return cwl.ArraySchema(items=new_type, type_="array")
    elif linkMerge == "merge_flattened":
        return merge_flatten_type(new_type)
    elif isinstance(sourcenames, list) and len(sourcenames) > 1:
        return cwl.ArraySchema(items=new_type, type_="array")
    else:
        return new_type


def param_for_source_id(
    process: Union[cwl.CommandLineTool, cwl.Workflow, cwl.ExpressionTool],
    sourcenames: Union[str, list[str]],
    parent: Optional[cwl.Workflow] = None,
    scatter_context: Optional[list[Optional[tuple[int, str]]]] = None,
) -> Union[list[cwl.InputParameter], cwl.InputParameter]:
    """Find the process input parameter that matches one of the given sourcenames."""
    if isinstance(sourcenames, str):
        sourcenames = [sourcenames]
    params: list[cwl.InputParameter] = []
    for sourcename in sourcenames:
        if not isinstance(process, cwl.Workflow):
            for param in process.inputs:
                if param.id.split("#")[-1] == sourcename.split("#")[-1]:
                    params.append(param)
                    if scatter_context is not None:
                        scatter_context.append(None)
        targets = [process]
        if parent:
            targets.append(parent)
        for target in targets:
            if isinstance(target, cwl.Workflow):
                for inp in target.inputs:
                    if inp.id.split("#")[-1] == sourcename.split("#")[-1]:
                        params.append(inp)
                        if scatter_context is not None:
                            scatter_context.append(None)
                for step in target.steps:
                    if (
                        "/".join(sourcename.split("#")[-1].split("/")[:-1])
                        == step.id.split("#")[-1]
                        and step.out
                    ):
                        step_run = cwl_utils.parser.utils.load_step(step)
                        cwl_utils.parser.utils.convert_stdstreams_to_files(step_run)
                        for outp in step.out:
                            outp_id = outp if isinstance(outp, str) else outp.id
                            if (
                                outp_id.split("#")[-1].split("/")[-1]
                                == sourcename.split("#")[-1].split("/")[-1]
                            ):
                                if step_run and step_run.outputs:
                                    for output in step_run.outputs:
                                        if (
                                            output.id.split("#")[-1].split("/")[-1]
                                            == sourcename.split("#")[-1].split("/")[-1]
                                        ):
                                            params.append(output)
                                            if scatter_context is not None:
                                                if isinstance(step.scatter, str):
                                                    scatter_context.append(
                                                        (
                                                            1,
                                                            step.scatterMethod
                                                            or "dotproduct",
                                                        )
                                                    )
                                                elif isinstance(
                                                    step.scatter, MutableSequence
                                                ):
                                                    scatter_context.append(
                                                        (
                                                            len(step.scatter),
                                                            step.scatterMethod
                                                            or "dotproduct",
                                                        )
                                                    )
                                                else:
                                                    scatter_context.append(None)
    if len(params) == 1:
        return params[0]
    elif len(params) > 1:
        return params
    raise WorkflowException(
        "param {} not found in {}\n{}.".format(
            sourcename,
            yaml_dumps(cwl.save(process)),
            (f" or\n {yaml_dumps(cwl.save(parent))}" if parent is not None else ""),
        )
    )
