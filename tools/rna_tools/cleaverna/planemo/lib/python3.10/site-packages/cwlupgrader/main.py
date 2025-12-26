#!/usr/bin/env python
"""Transforms draft-3 CWL documents into v1.0 as idiomatically as possible."""

import argparse
import copy
import logging
import os
import os.path
import stat
import sys
from collections.abc import MutableMapping, MutableSequence, Sequence
from pathlib import Path
from typing import Any, Callable, Optional, Union

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap  # for consistent sort order
from schema_salad.sourceline import SourceLine, add_lc_filename, cmap

_logger = logging.getLogger("cwl-upgrader")  # pylint: disable=invalid-name
defaultStreamHandler = logging.StreamHandler()  # pylint: disable=invalid-name
_logger.addHandler(defaultStreamHandler)
_logger.setLevel(logging.INFO)

yaml = ruamel.yaml.main.YAML(typ="rt")
yaml.allow_duplicate_keys = True
yaml.preserve_quotes = True
yaml.default_flow_style = False


def parse_args(args: list[str]) -> argparse.Namespace:
    """Argument parser."""
    parser = argparse.ArgumentParser(
        description="Tool to upgrade CWL documents from one version to another. "
        "Supports upgrading 'draft-3', 'v1.0', and 'v1.1' to 'v1.2'",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--v1-only", help="Don't upgrade past cwlVersion: v1.0", action="store_true"
    )
    parser.add_argument(
        "--v1.1-only",
        dest="v1_1_only",
        help="Don't upgrade past cwlVersion: v1.1",
        action="store_true",
    )
    parser.add_argument(
        "--dir", help="Directory in which to save converted files", default=Path.cwd()
    )
    parser.add_argument(
        "--always-write",
        help="Always write a file, even if no changes were made.",
        action="store_true",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more CWL documents.",
    )
    return parser.parse_args(args)


def main(args: Optional[list[str]] = None) -> int:
    """Hook to set the args."""
    if not args:
        args = sys.argv[1:]
    return run(parse_args(args))


def run(args: argparse.Namespace) -> int:
    """Main function."""
    imports: set[str] = set()
    if args.dir and not os.path.exists(args.dir):
        os.makedirs(args.dir)
    for path in args.inputs:
        _logger.info("Processing %s", path)
        document = load_cwl_document(path)
        if "cwlVersion" not in document:
            _logger.warn("No cwlVersion found in %s, skipping it.", path)
        else:
            if document["cwlVersion"] == "v1.0":
                if args.v1_only:
                    _logger.info("Skipping v1.0 document as requested: %s.", path)
                    continue
            elif document["cwlVersion"] == "v1.1":
                if args.v1_1_only:
                    _logger.info("Skipping v1.1 document as requested: %s.", path)
                    continue

            if args.v1_only:
                target_version = "v1.0"
            elif args.v1_1_only:
                target_version = "v1.1"
            else:
                target_version = "latest"
            upgraded_document = upgrade_document(
                document,
                args.dir,
                target_version=target_version,
                imports=imports,
            )
            if upgraded_document is not document or not args.always_write:
                write_cwl_document(upgraded_document, Path(path).name, args.dir)
    return 0


def upgrade_document(
    document: Any,
    output_dir: str,
    target_version: Optional[str] = "latest",
    imports: Optional[set[str]] = None,
) -> Any:
    if imports is None:
        imports = set()
    supported_versions = ["v1.0", "v1.1", "v1.2", "latest"]
    if target_version not in supported_versions:
        _logger.error(f"Unsupported target cwlVersion: {target_version}")
        return

    version = document["cwlVersion"]
    main_updater = None
    inner_updater = None

    if version == "cwl:draft-3" or version == "draft-3":
        if target_version == "v1.0":
            main_updater = draft3_to_v1_0
            inner_updater = _draft3_to_v1_0
        elif target_version == "v1.1":
            main_updater = draft3_to_v1_1
            inner_updater = _draft3_to_v1_1
        elif target_version == "v1.2":
            main_updater = draft3_to_v1_2
            inner_updater = _draft3_to_v1_2
        elif target_version == "latest":
            main_updater = draft3_to_v1_2
            inner_updater = _draft3_to_v1_2
    elif version == "v1.0":
        if target_version == "v1.0":
            _logger.info("Not upgrading v1.0 document as requested.")
            return
        elif target_version == "v1.1":
            main_updater = v1_0_to_v1_1
            inner_updater = _v1_0_to_v1_1
        elif target_version == "v1.2":
            main_updater = v1_0_to_v1_2
            inner_updater = _v1_0_to_v1_2
        elif target_version == "latest":
            main_updater = v1_0_to_v1_2
            inner_updater = _v1_0_to_v1_2
    elif version == "v1.1":
        if target_version == "v1.1":
            _logger.info("Not upgrading v1.1 document as requested.")
            return
        elif target_version == "v1.2":
            main_updater = v1_1_to_v1_2
            inner_updater = _v1_1_to_v1_2
        elif target_version == "latest":
            main_updater = v1_1_to_v1_2
            inner_updater = _v1_1_to_v1_2
    elif version == "v1.2":
        if target_version == "v1.2":
            _logger.info("Not upgrading v1.2 document as requested.")
            return document
        elif target_version == "latest":
            return document
    else:
        _logger.error(f"Unknown cwlVersion in source document: {version}")
        return

    if main_updater is None or inner_updater is None:
        _logger.error(f"Cannot downgrade from cwlVersion {version} to {target_version}")
        return

    process_imports(document, imports, inner_updater, output_dir)
    return main_updater(document, output_dir)


def load_cwl_document(path: str) -> Any:
    """
    Load the given path using the Ruamel YAML round-trip loader.

    Also ensures that the filename is recorded so that SourceLine can produce
    informative error messages.
    """
    with open(path) as entry:
        document = yaml.load(entry)
        add_lc_filename(document, entry.name)
    return document


def write_cwl_document(document: Any, name: str, dirname: str) -> None:
    """
    Serialize the document using the Ruamel YAML round trip dumper.

    Will also prepend "#!/usr/bin/env cwl-runner\n" and
    set the executable bit if it is a CWL document.
    """
    ruamel.yaml.scalarstring.walk_tree(document)
    path = Path(dirname) / name
    with open(path, "w") as handle:
        if "cwlVersion" in document:
            if not (
                document.ca
                and document.ca.comment
                and "cwl-runner" in document.ca.comment[1][0].value
            ):
                handle.write("#!/usr/bin/env cwl-runner\n")
        yaml.dump(document, stream=handle)
    if "cwlVersion" in document:
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def process_imports(
    document: Any, imports: set[str], updater: Callable[[Any, str], Any], outdir: str
) -> None:
    """Find any '$import's and process them."""
    if isinstance(document, CommentedMap):
        for key, value in document.items():
            if key == "$import":
                if value not in imports:
                    write_cwl_document(
                        updater(
                            load_cwl_document(
                                Path(document.lc.filename).parent / value
                            ),
                            outdir,
                        ),
                        Path(value).name,
                        outdir,
                    )
                    imports.add(value)
            else:
                process_imports(value, imports, updater, outdir)
    elif isinstance(document, MutableSequence):
        for entry in document:
            process_imports(entry, imports, updater, outdir)


def v1_0_to_v1_1(document: CommentedMap, outdir: str) -> CommentedMap:
    """CWL v1.0.x to v1.1 transformation loop."""
    _v1_0_to_v1_1(document, outdir)
    for key, value in document.items():
        with SourceLine(document, key, Exception):
            if isinstance(value, CommentedMap):
                document[key] = _v1_0_to_v1_1(value, outdir)
            elif isinstance(value, list):
                for index, entry in enumerate(value):
                    if isinstance(entry, CommentedMap):
                        value[index] = _v1_0_to_v1_1(entry, outdir)
    document["cwlVersion"] = "v1.1"
    return sort_v1_0(document)


def v1_0_to_v1_2(document: CommentedMap, outdir: str) -> CommentedMap:
    """CWL v1.0.x to v1.2 transformation."""
    document = v1_0_to_v1_1(document, outdir)
    document = v1_1_to_v1_2(document, outdir)
    return document


def v1_1_to_v1_2(document: CommentedMap, outdir: str) -> CommentedMap:
    """CWL v1.1 to v1.2 transformation."""
    document = _v1_1_to_v1_2(document, outdir)
    document["cwlVersion"] = "v1.2"
    return document


def draft3_to_v1_0(document: CommentedMap, outdir: str) -> CommentedMap:
    """Transformation loop."""
    _draft3_to_v1_0(document, outdir)
    if isinstance(document, MutableMapping):
        for key, value in document.items():
            with SourceLine(document, key, Exception):
                if isinstance(value, CommentedMap):
                    document[key] = _draft3_to_v1_0(value, outdir)
                elif isinstance(value, list):
                    for index, entry in enumerate(value):
                        if isinstance(entry, CommentedMap):
                            value[index] = _draft3_to_v1_0(entry, outdir)
    document["cwlVersion"] = "v1.0"
    return sort_v1_0(document)


def draft3_to_v1_1(document: CommentedMap, outdir: str) -> CommentedMap:
    """transformation loop."""
    return v1_0_to_v1_1(draft3_to_v1_0(document, outdir), outdir)


def draft3_to_v1_2(document: CommentedMap, outdir: str) -> CommentedMap:
    """transformation loop."""
    return v1_1_to_v1_2(v1_0_to_v1_1(draft3_to_v1_0(document, outdir), outdir), outdir)


def _draft3_to_v1_0(document: CommentedMap, outdir: str) -> CommentedMap:
    """Inner loop for transforming draft-3 to v1.0."""
    if "class" in document:
        if document["class"] == "Workflow":
            workflow_clean(document)
        elif document["class"] == "File":
            document["location"] = document.pop("path")
        elif document["class"] == "CommandLineTool":
            input_output_clean(document)
            hints_and_requirements_clean(document)
            if (
                isinstance(document["baseCommand"], list)
                and len(document["baseCommand"]) == 1
            ):
                document["baseCommand"] = document["baseCommand"][0]
            if "arguments" in document and not document["arguments"]:
                del document["arguments"]
    clean_secondary_files(document)

    if "description" in document:
        document["doc"] = document.pop("description")

    return document


def _draft3_to_v1_1(document: CommentedMap, outdir: str) -> CommentedMap:
    return v1_0_to_v1_1(_draft3_to_v1_0(document, outdir), outdir)


def _draft3_to_v1_2(document: CommentedMap, outdir: str) -> CommentedMap:
    return _draft3_to_v1_1(document, outdir)  # nothing needs doing for 1.2


WORKFLOW_INPUT_INPUTBINDING = (
    "{}[cwl-upgrader_v1_0_to_v1_1] Original input had the following "
    "(unused) inputBinding element: {}"
)

V1_0_TO_V1_1_REWRITE = {
    "http://commonwl.org/cwltool#WorkReuse": "WorkReuse",
    "http://arvados.org/cwl#ReuseRequirement": "WorkReuse",
    "http://commonwl.org/cwltool#TimeLimit": "ToolTimeLimit",
    "http://commonwl.org/cwltool#NetworkAccess": "NetworkAccess",
    "http://commonwl.org/cwltool#InplaceUpdateRequirement": "InplaceUpdateRequirement",
    "http://commonwl.org/cwltool#LoadListingRequirement": "LoadListingRequirement",
}


def _v1_0_to_v1_1(document: CommentedMap, outdir: str) -> CommentedMap:
    """Inner loop for transforming draft-3 to v1.0."""
    if "class" in document:
        if document["class"] == "Workflow":
            upgrade_v1_0_hints_and_reqs(document)
            move_up_loadcontents(document)
            cleanup_v1_0_input_bindings(document)
            steps = document["steps"]
            if isinstance(steps, MutableSequence):
                for index, entry in enumerate(steps):
                    with SourceLine(steps, index, Exception):
                        upgrade_v1_0_hints_and_reqs(entry)
                        if "run" in entry and isinstance(entry["run"], CommentedMap):
                            process = entry["run"]
                            _v1_0_to_v1_1(process, outdir)
                            if "cwlVersion" in process:
                                del process["cwlVersion"]
                        elif isinstance(entry["run"], str) and "#" not in entry["run"]:
                            path = Path(document.lc.filename).parent / entry["run"]
                            process = v1_0_to_v1_1(load_cwl_document(str(path)), outdir)
                            write_cwl_document(process, path.name, outdir)
            elif isinstance(steps, MutableMapping):
                for step_name in steps:
                    with SourceLine(steps, step_name, Exception):
                        entry = steps[step_name]
                        upgrade_v1_0_hints_and_reqs(entry)
                        if "run" in entry:
                            if isinstance(entry["run"], CommentedMap):
                                process = entry["run"]
                                _v1_0_to_v1_1(process, outdir)
                                if "cwlVersion" in process:
                                    del process["cwlVersion"]
                            elif (
                                isinstance(entry["run"], str)
                                and "#" not in entry["run"]
                            ):
                                path = Path(document.lc.filename).parent / entry["run"]
                                process = v1_0_to_v1_1(
                                    load_cwl_document(str(path)), outdir
                                )
                                write_cwl_document(process, path.name, outdir)
                            elif isinstance(entry["run"], str) and "#" in entry["run"]:
                                pass  # reference to $graph entry
                            else:
                                raise Exception(
                                    "'run' entry was neither a CWL Process nor "
                                    "a path to one: %s.",
                                    entry["run"],
                                )
        elif document["class"] == "CommandLineTool":
            upgrade_v1_0_hints_and_reqs(document)
            move_up_loadcontents(document)
            network_access = has_hint_or_req(document, "NetworkAccess")
            listing = has_hint_or_req(document, "LoadListingRequirement")
            reqs = document.get("requirements", {})
            # TODO: add comments to explain the extra hints
            if isinstance(reqs, MutableSequence):
                if not network_access:
                    reqs.append({"class": "NetworkAccess", "networkAccess": True})
                if not listing:
                    reqs.append(
                        cmap(
                            {
                                "class": "LoadListingRequirement",
                                "loadListing": "deep_listing",
                            }
                        )
                    )
            elif isinstance(reqs, MutableMapping):
                if not network_access:
                    reqs["NetworkAccess"] = {"networkAccess": True}
                if not listing:
                    reqs["LoadListingRequirement"] = cmap(
                        {"loadListing": "deep_listing"}
                    )
            if "requirements" not in document:
                document["requirements"] = reqs
        elif document["class"] == "ExpressionTool":
            move_up_loadcontents(document)
            cleanup_v1_0_input_bindings(document)
    return document


def _v1_0_to_v1_2(document: CommentedMap, outdir: str) -> CommentedMap:
    document = _v1_0_to_v1_1(document, outdir)
    return _v1_1_to_v1_2(document, outdir)


def _v1_1_to_v1_2(document: CommentedMap, outdir: str) -> CommentedMap:
    if "class" in document:
        if document["class"] == "Workflow":
            steps = document["steps"]
            if isinstance(steps, MutableSequence):
                for index, entry in enumerate(steps):
                    with SourceLine(steps, index, Exception):
                        if "run" in entry and isinstance(entry["run"], CommentedMap):
                            process = entry["run"]
                            _v1_1_to_v1_2(process, outdir)
                            if "cwlVersion" in process:
                                del process["cwlVersion"]

                        elif isinstance(entry["run"], str) and "#" not in entry["run"]:
                            if hasattr(document.lc, "filename"):
                                dirname = Path(document.lc.filename).parent
                            else:
                                dirname = Path(outdir)
                            path = dirname / entry["run"]
                            process = v1_1_to_v1_2(load_cwl_document(str(path)), outdir)
                            write_cwl_document(process, path.name, outdir)
            elif isinstance(steps, MutableMapping):
                for step_name in steps:
                    with SourceLine(steps, step_name, Exception):
                        entry = steps[step_name]
                        if "run" in entry:
                            if isinstance(entry["run"], CommentedMap):
                                process = entry["run"]
                                _v1_1_to_v1_2(process, outdir)
                                if "cwlVersion" in process:
                                    del process["cwlVersion"]
                            elif (
                                isinstance(entry["run"], str)
                                and "#" not in entry["run"]
                            ):
                                if hasattr(document.lc, "filename"):
                                    dirname = Path(document.lc.filename).parent
                                else:
                                    dirname = Path(outdir)
                                path = dirname / entry["run"]
                                process = v1_1_to_v1_2(
                                    load_cwl_document(str(path)), outdir
                                )
                                write_cwl_document(process, path.name, outdir)
                            elif isinstance(entry["run"], str) and "#" in entry["run"]:
                                pass  # reference to $graph entry
                            else:
                                raise Exception(
                                    "'run' entry was neither a CWL Process nor "
                                    "a path to one: %s.",
                                    entry["run"],
                                )
    return document


def cleanup_v1_0_input_bindings(document: dict[str, Any]) -> None:
    """In v1.1 Workflow or ExpressionTool level inputBindings are deprecated."""

    def cleanup(inp: dict[str, Any]) -> None:
        """Serialize non loadContents fields and add that to the doc."""
        if "inputBinding" in inp:
            bindings = inp["inputBinding"]
            for field in list(bindings.keys()):
                if field != "loadContents":
                    prefix = "" if "doc" not in inp else "{}\n".format(inp["doc"])
                    inp["doc"] = WORKFLOW_INPUT_INPUTBINDING.format(prefix, field)
                    del bindings[field]
            if not bindings:
                del inp["inputBinding"]

    inputs = document["inputs"]
    if isinstance(inputs, MutableSequence):
        for entry in inputs:
            cleanup(entry)
    elif isinstance(inputs, MutableMapping):
        for input_name in inputs:
            cleanup(inputs[input_name])


def move_up_loadcontents(document: dict[str, Any]) -> None:
    """Promote 'loadContents' up a level for CWL v1.1."""

    def cleanup(inp: dict[str, Any]) -> None:
        """Move loadContents to the preferred location."""
        if "inputBinding" in inp:
            bindings = inp["inputBinding"]
            for field in list(bindings.keys()):
                if field == "loadContents":
                    inp[field] = bindings.pop(field)

    inputs = document["inputs"]
    if isinstance(inputs, MutableSequence):
        for entry in inputs:
            cleanup(entry)
    elif isinstance(inputs, MutableMapping):
        for input_name in inputs:
            cleanup(inputs[input_name])


def upgrade_v1_0_hints_and_reqs(document: dict[str, Any]) -> None:
    """Rename some pre-v1.1 extensions to their official CWL v1.1 names."""
    for extra in ("requirements", "hints"):
        if extra in document:
            with SourceLine(document, extra, Exception):
                if isinstance(document[extra], MutableMapping):
                    for req_name in document[extra]:
                        with SourceLine(document[extra], req_name, Exception):
                            if req_name in V1_0_TO_V1_1_REWRITE:
                                document[extra][V1_0_TO_V1_1_REWRITE[req_name]] = (
                                    document[extra].pop(req_name)
                                )
                elif isinstance(document[extra], MutableSequence):
                    for index, entry in enumerate(document[extra]):
                        with SourceLine(document[extra], index, Exception):
                            if (
                                isinstance(entry, MutableMapping)
                                and "class" in entry
                                and entry["class"] in V1_0_TO_V1_1_REWRITE
                            ):
                                entry["class"] = V1_0_TO_V1_1_REWRITE[entry["id"]]
                else:
                    raise Exception(
                        "{} section must be either a list of dictionaries "
                        "or a dictionary of dictionaries!: {}".format(
                            extra, document[extra]
                        )
                    )


def has_hint_or_req(document: dict[str, Any], name: str) -> bool:
    """Detects an existing named hint or requirement."""
    for extra in ("requirements", "hints"):
        if extra in document:
            with SourceLine(document, extra, Exception):
                if isinstance(document[extra], MutableMapping):
                    if name in document[extra]:
                        return True
                elif isinstance(document[extra], MutableSequence):
                    for index, entry in enumerate(document[extra]):
                        with SourceLine(document[extra], index, Exception):
                            if "class" == entry and entry["class"] == name:
                                return True
    return False


def workflow_clean(document: dict[str, Any]) -> None:
    """Transform draft-3 style Workflows to more idiomatic v1.0"""
    input_output_clean(document)
    hints_and_requirements_clean(document)
    outputs = document["outputs"]
    for index, output_id in enumerate(outputs):
        with SourceLine(outputs, index, Exception):
            outputs[output_id]["outputSource"] = (
                outputs[output_id].pop("source").lstrip("#").replace(".", "/")
            )
    new_steps = CommentedMap()
    for index, step in enumerate(document["steps"]):
        with SourceLine(document["steps"], index, Exception):
            new_step = CommentedMap()
            new_step.update(step)
            step = new_step
            step_id = step.pop("id")
            step_id_len = len(step_id) + 1
            step["out"] = []
            for index2, outp in enumerate(step["outputs"]):
                with SourceLine(step["outputs"], index2, Exception):
                    clean_outp_id = outp["id"]
                    if clean_outp_id.startswith(step_id):
                        clean_outp_id = clean_outp_id[step_id_len:]
                    step["out"].append(clean_outp_id)
            del step["outputs"]
            ins = CommentedMap()
            for index3, inp in enumerate(step["inputs"]):
                with SourceLine(step["inputs"], index3, Exception):
                    ident = inp["id"]
                    if ident.startswith(step_id):
                        ident = ident[step_id_len:]
                    if "source" in inp:
                        with SourceLine(inp, "source", Exception):
                            if isinstance(inp["source"], str):
                                inp["source"] = (
                                    inp["source"].lstrip("#").replace(".", "/")
                                )
                            else:
                                for index4, inp_source in enumerate(inp["source"]):
                                    with SourceLine(inp["source"], index4, Exception):
                                        inp["source"][index4] = inp_source.lstrip(
                                            "#"
                                        ).replace(".", "/")
                    del inp["id"]
                    if len(inp) > 1:
                        ins[ident] = inp
                    elif len(inp) == 1:
                        if "source" in inp:
                            ins[ident] = inp.popitem()[1]
                        else:
                            ins[ident] = inp
                    else:
                        ins[ident] = {}
            step["in"] = ins
            del step["inputs"]
            if "scatter" in step:
                with SourceLine(step, "scatter", Exception):
                    if isinstance(step["scatter"], str) == 1:
                        source = step["scatter"]
                        if source.startswith(step_id):
                            source = source[step_id_len:]
                        step["scatter"] = source
                    elif isinstance(step["scatter"], list) and len(step["scatter"]) > 1:
                        step["scatter"] = []
                        for index4, source in enumerate(step["scatter"]):
                            with SourceLine(step["scatter"], index4, Exception):
                                if source.startswith(step_id):
                                    source = source[step_id_len:]
                                step["scatter"].append(source)
                    else:
                        source = step["scatter"][0]
                        if source.startswith(step_id):
                            source = source[step_id_len:]
                        step["scatter"] = source
            if "description" in step:
                step["doc"] = step.pop("description")
            new_steps[step_id.lstrip("#")] = step
    document["steps"] = new_steps


def input_output_clean(document: dict[str, Any]) -> None:
    """Transform draft-3 style input/output listings into idiomatic v1.0."""
    for param_type in ["inputs", "outputs"]:
        if param_type not in document:
            break
        new_section = CommentedMap()
        meta = False
        for index, param in enumerate(document[param_type]):
            with SourceLine(document[param_type], index, Exception):
                if "$import" in param:
                    meta = True
        if not meta:
            for index2, param2 in enumerate(document[param_type]):
                with SourceLine(document[param_type], index2, Exception):
                    param_id = param2.pop("id").lstrip("#")
                    if "type" in param2:
                        param2["type"] = shorten_type(param2["type"])
                        array_type_raise_sf(param2)
                    if "description" in param2:
                        param2["doc"] = param2.pop("description")
                    if len(param2) > 1:
                        new_section[param_id] = sort_input_or_output(param2)
                    elif "type" in param2 and isinstance(param2["type"], str):
                        new_section[param_id] = param2.popitem()[1]
                    else:
                        new_section[param_id] = param2
            document[param_type] = new_section


def array_type_raise_sf(param: MutableMapping[str, Any]) -> None:
    """Move up draft-3 secondaryFile specs on File members in Arrays."""
    typ = param["type"]
    if isinstance(typ, MutableSequence):
        for index, param2 in enumerate(typ):
            with SourceLine(typ, index, Exception):
                if isinstance(param2, MutableMapping) and "type" in param2:
                    array_type_raise_sf(param2)
    elif (
        isinstance(typ, MutableMapping)
        and "type" in typ
        and typ["type"] == "array"
        and "items" in typ
        and "File" in typ["items"]
        and "secondaryFiles" in typ
    ):
        param["secondaryFiles"] = typ["secondaryFiles"]
        del typ["secondaryFiles"]


def hints_and_requirements_clean(document: dict[str, Any]) -> None:
    """Transform draft-3 style hints/reqs into idiomatic v1.0 hints/reqs."""
    for section in ["hints", "requirements"]:
        if section in document:
            new_section = {}
            meta = False
            for index, entry in enumerate(document[section]):
                with SourceLine(document[section], index, Exception):
                    if isinstance(entry, MutableMapping):
                        if "$import" in entry or "$include" in entry:
                            meta = True
            for index2, entry2 in enumerate(document[section]):
                with SourceLine(document[section], index2, Exception):
                    if isinstance(entry2, MutableMapping):
                        if (
                            "class" in entry2
                            and entry2["class"] == "CreateFileRequirement"
                        ):
                            entry2["class"] = "InitialWorkDirRequirement"
                            entry2["listing"] = []
                            for filedef in entry2["fileDef"]:
                                entry2["listing"].append(
                                    {
                                        "entryname": filedef["filename"],
                                        "entry": filedef["fileContent"],
                                    }
                                )
                            del entry2["fileDef"]
                    if not meta:
                        new_section[entry2["class"]] = entry2
                        del entry2["class"]
            if not meta:
                document[section] = new_section


def shorten_type(type_obj: Union[str, list[Any]]) -> Union[str, list[Any]]:
    """Transform draft-3 style type declarations into idiomatic v1.0 types."""
    if isinstance(type_obj, str) or not isinstance(type_obj, Sequence):
        return type_obj
    new_type: list[str] = []
    for entry in type_obj:  # find arrays that we can shorten and do so
        if isinstance(entry, dict):
            if entry["type"] == "array" and isinstance(entry["items"], str):
                entry = entry["items"] + "[]"
            elif entry["type"] == "enum":
                entry = sort_enum(entry)
        new_type.extend([entry])
    if len(new_type) == 2:
        if "null" in new_type:
            type_copy = copy.deepcopy(new_type)
            type_copy.remove("null")
            if isinstance(type_copy[0], str):
                return type_copy[0] + "?"
    if len(new_type) == 1:
        return new_type[0]
    return new_type


def clean_secondary_files(document: dict[str, Any]) -> None:
    """Cleanup for secondaryFiles"""
    if "secondaryFiles" in document:
        for i, sfile in enumerate(document["secondaryFiles"]):
            if "$(" in sfile or "${" in sfile:
                document["secondaryFiles"][i] = sfile.replace(
                    '"path"', '"location"'
                ).replace(".path", ".location")


def sort_v1_0(document: dict[str, Any]) -> CommentedMap:
    """Sort the sections of the CWL document in a more meaningful order."""
    keyorder = [
        "cwlVersion",
        "class",
        "id",
        "label",
        "doc",
        "requirements",
        "hints",
        "inputs",
        "stdin",
        "baseCommand",
        "steps",
        "expression",
        "arguments",
        "stderr",
        "stdout",
        "outputs",
        "successCodes",
        "temporaryFailCodes",
        "permanentFailCodes",
    ]
    return CommentedMap(
        sorted(
            document.items(),
            key=lambda i: keyorder.index(i[0]) if i[0] in keyorder else 100,
        )
    )


def sort_enum(enum: dict[str, Any]) -> dict[str, Any]:
    """Sort the enum type definitions in a more meaningful order."""
    keyorder = ["type", "name", "label", "symbols", "inputBinding"]
    return CommentedMap(
        sorted(
            enum.items(),
            key=lambda i: keyorder.index(i[0]) if i[0] in keyorder else 100,
        )
    )


def sort_input_or_output(io_def: dict[str, Any]) -> dict[str, Any]:
    """Sort the input definitions in a more meaningful order."""
    keyorder = [
        "label",
        "doc",
        "type",
        "format",
        "secondaryFiles",
        "default",
        "inputBinding",
        "outputBinding",
        "streamable",
    ]
    return CommentedMap(
        sorted(
            io_def.items(),
            key=lambda i: keyorder.index(i[0]) if i[0] in keyorder else 100,
        )
    )


if __name__ == "__main__":
    sys.exit(main())
