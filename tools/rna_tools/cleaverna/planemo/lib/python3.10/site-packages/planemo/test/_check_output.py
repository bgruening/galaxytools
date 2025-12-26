"""Check an output file from a generalize artifact test."""

import os
import tempfile

import requests
from galaxy.tool_util.parser.interface import TestCollectionOutputDef
from galaxy.tool_util.verify import verify
from galaxy.tool_util.verify.interactor import verify_collection
from galaxy.util import unicodify


def check_output(runnable, output_properties, test_properties, **kwds):
    """Use galaxy-tool-util to check a test output.

    Return a list of strings describing the problems encountered,
    and empty list indicates no problems were detected.

    Currently this will only ever return at most one detected problem because
    of the way galaxy-tool-util throws exceptions instead of returning individual
    descriptions - but this may be enhanced in the future.
    """
    checker = _check_output_collection if for_collections(test_properties) else _check_output_file
    return checker(runnable, output_properties, test_properties, **kwds)


def for_collections(test_properties):
    return "element_tests" in test_properties


def _check_output_collection(runnable, output_properties, test_properties, **kwds):
    data_collection = output_properties

    output_def = TestCollectionOutputDef.from_dict(test_properties)

    def verify_dataset(element, element_attrib, element_outfile):
        if element_outfile:
            element_attrib["path"] = element_outfile
        _verify_output_file(runnable, element["_output_object"], element_attrib)

    problems = []
    try:
        verify_collection(output_def, data_collection, verify_dataset)
    except AssertionError as e:
        problems.append(unicodify(e))

    return problems


def _verify_output_file(runnable, output_properties, test_properties, **kwds):
    get_filename = _test_filename_getter(runnable)
    path = output_properties["path"]
    with open(path, "rb") as fh:
        output_content = fh.read()
    # Support Galaxy-like file location (using "file") or CWL-like ("path" or "location").
    expected_file = test_properties.get("file", None)
    if expected_file is None:
        expected_file = test_properties.get("path", None)
    if expected_file is None:
        location = test_properties.get("location")
        if location:
            if location.startswith(("http://", "https://")):
                expected_file = _get_location(location)
            else:
                expected_file = location.split("file://", 1)[-1]

    test_data_target_dir = kwds.get("test_data_target_dir", None)
    item_label = "Output with path %s" % path
    if "asserts" in test_properties:
        # TODO: break fewer abstractions here...
        from galaxy.tool_util.parser.yaml import __to_test_assert_list

        test_properties["assert_list"] = __to_test_assert_list(test_properties["asserts"])
    verify(
        item_label,
        output_content,
        attributes=test_properties,
        filename=expected_file,
        get_filename=get_filename,
        keep_outputs_dir=test_data_target_dir,
        verify_extra_files=None,
    )


def _check_output_file(runnable, output_properties, test_properties, **kwds):
    problems = []
    try:
        _verify_output_file(runnable, output_properties, test_properties, **kwds)
    except AssertionError as e:
        problems.append(unicodify(e))

    return problems


def _get_location(location):
    data_file = tempfile.NamedTemporaryFile(prefix="planemo_test_file_", delete=False)
    with requests.get(location, stream=True) as r:
        r.raise_for_status()

        for chunk in r.iter_content():
            if chunk:
                data_file.write(chunk)
        return data_file.name


def _test_filename_getter(runnable):
    def get_filename(name):
        artifact_directory = os.path.dirname(runnable.path)
        return os.path.join(artifact_directory, name)

    return get_filename


__all__ = ("check_output",)
