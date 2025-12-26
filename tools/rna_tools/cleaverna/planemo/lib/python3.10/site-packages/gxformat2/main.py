"""Module containing :func:`convert_and_import_workflow`."""
import os

from .converter import python_to_workflow, yaml_to_workflow
from .interface import BioBlendImporterGalaxyInterface
from .yaml import ordered_load


def convert_and_import_workflow(has_workflow, **kwds):
    """Conversion and import Format 2 workflows into a target Galaxy instance."""
    galaxy_interface = kwds.get("galaxy_interface", None)
    if galaxy_interface is None:
        galaxy_interface = BioBlendImporterGalaxyInterface(**kwds)

    source_type = kwds.get("source_type", None)
    workflow_directory = kwds.get("workflow_directory", None)
    if source_type == "path":
        workflow_path = has_workflow
        if workflow_directory is None:
            workflow_directory = os.path.dirname(has_workflow)
        with open(workflow_path) as f:
            has_workflow = ordered_load(f)

    if workflow_directory is not None:
        workflow_directory = os.path.abspath(workflow_directory)

    convert = kwds.get("convert", True)
    raw_yaml = kwds.get("raw_yaml", False)
    if raw_yaml and convert:
        raise Exception("Incompatible options selected.")
    if convert:
        if isinstance(has_workflow, dict):
            workflow = python_to_workflow(has_workflow, galaxy_interface, workflow_directory)
        else:
            workflow = yaml_to_workflow(has_workflow, galaxy_interface, workflow_directory)
    else:
        workflow = has_workflow
        if not isinstance(workflow, dict) and not raw_yaml:
            workflow = ordered_load(workflow)
        else:
            workflow = {"yaml_content": workflow}

    name = kwds.get("name", None)
    if name is not None:
        workflow["name"] = name
    publish = kwds.get("publish", False)
    exact_tools = kwds.get("exact_tools", False)
    fill_defaults = kwds.get("fill_defaults", True)
    import_kwds = {
        "fill_defaults": fill_defaults
    }
    if publish:
        import_kwds["publish"] = True
    if exact_tools:
        import_kwds["exact_tools"] = True
    return galaxy_interface.import_workflow(workflow, **import_kwds)


__all__ = (
    'convert_and_import_workflow',
)
