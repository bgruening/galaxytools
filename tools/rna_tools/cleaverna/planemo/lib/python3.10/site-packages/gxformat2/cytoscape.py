"""Build standalone visualization for Galaxy workflows."""
import argparse
import json
import os
import string
import sys

from gxformat2.model import ensure_step_position
from gxformat2.normalize import steps_normalized

CYTOSCAPE_JS_TEMPLATE = os.path.join(os.path.dirname(__file__), 'cytoscape.html')
MAIN_TS_PREFIX = "toolshed.g2.bx.psu.edu/repos/"
SCRIPT_DESCRIPTION = """
This script converts an executable Galaxy workflow (in either format - Format 2
or native .ga) into a format for visualization with Cytoscape
(https://cytoscape.org/).

If the target output path ends with .html this script will output a HTML
page with the workflow visualized using cytoscape.js. Otherwise, this script
will output a JSON description of the workflow elements for consumption by
Cytoscape.
"""


def to_cytoscape(workflow_path: str, output_path=None):
    """Produce cytoscape output for supplied workflow path."""
    if output_path is None:
        output_path, _ = os.path.splitext(workflow_path)
        output_path += ".html"

    steps = steps_normalized(workflow_path=workflow_path)
    elements = []
    for i, step in enumerate(steps):
        step_id = step.get("id") or step.get("label") or str(i)
        step_type = step.get("type") or 'tool'
        classes = [f"type_{step_type}"]
        if step_type in ['tool', 'subworkflow']:
            classes.append("runnable")
        else:
            classes.append("input")

        tool_id = step.get("tool_id")
        if tool_id and tool_id.startswith(MAIN_TS_PREFIX):
            tool_id = tool_id[len(MAIN_TS_PREFIX):]
        label = step.get("id") or step.get("label") or (f"tool:{tool_id}" if tool_id else str(i))
        ensure_step_position(step, i)
        node_position = dict(x=int(step["position"]["left"]), y=int(step["position"]["top"]))
        repo_link = None
        if "tool_shed_repository" in step:
            repo = step["tool_shed_repository"]
            repo_link = "https://" + repo["tool_shed"] + "/view/" + repo["owner"] + "/" + repo["name"] + "/" + repo["changeset_revision"]
        node_data = {
            "id": step_id,
            "label": label,
            "doc": step.get("doc"),
            "tool_id": step.get("tool_id"),
            "step_type": step_type,
            "repo_link": repo_link
        }
        elements.append({"group": "nodes", "data": node_data, "classes": classes, "position": node_position})
        for key, value in (step.get("in") or {}).items():
            # handle lists?
            if isinstance(value, dict) and 'source' in value:
                value = value["source"]
            elif isinstance(value, dict):
                continue
            if "/" in value:
                from_step, output = value.split("/", 1)
            else:
                from_step, output = value, None
            edge_id = f"{step_id}__to__{from_step}"
            edge_data = {"id": edge_id, "source": from_step, "target": step_id, "input": key, "output": output}
            elements.append({"group": "edges", "data": edge_data})

    if output_path.endswith(".html"):
        with open(CYTOSCAPE_JS_TEMPLATE) as f:
            template = f.read()
        viz = string.Template(template).safe_substitute(elements=json.dumps(elements))
        with open(output_path, "w") as f:
            f.write(viz)
    else:
        with open(output_path, "w") as f:
            json.dump(elements, f)


def main(argv=None):
    """Entry point for building Cytoscape visualizations of Galaxy workflows."""
    if argv is None:
        argv = sys.argv[1:]

    args = _parser().parse_args(argv)

    workflow_path = args.input_path
    output_path = args.output_path
    to_cytoscape(workflow_path, output_path)


def _parser():
    parser = argparse.ArgumentParser(description=SCRIPT_DESCRIPTION)
    parser.add_argument('input_path', metavar='INPUT', type=str,
                        help='input workflow path (.ga/gxwf.yml)')
    parser.add_argument('output_path', metavar='OUTPUT', type=str, nargs="?",
                        help='output viz path (.json/.html)')
    return parser


if __name__ == "__main__":
    main()
