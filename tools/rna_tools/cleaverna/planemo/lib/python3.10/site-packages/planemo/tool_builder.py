"""This module contains :func:`build` to build tool descriptions.

This class is used by the `tool_init` command and can be used to build
Galaxy and CWL tool descriptions.
"""

import os
import re
import shlex
import shutil
import subprocess
from collections import namedtuple

from planemo import (
    io,
    templates,
)
from planemo.autopygen.argument_parser_conversion import (
    command_from_decoy,
    obtain_and_convert_parser,
    xml_from_decoy,
    xml_to_string,
)

REUSING_MACROS_MESSAGE = (
    "Macros file macros.xml already exists, assuming  it has relevant planemo-generated definitions."
)
DEFAULT_CWL_VERSION = "v1.0"

TOOL_TEMPLATE = """<tool id="{{id}}" name="{{name}}" version="{{version}}+galaxy0" python_template_version="3.5" profile="21.05">
{%- if description %}
    <description>{{ description }}</description>
{%- endif %}
{%- if macros %}
    <macros>
        <import>macros.xml</import>
    </macros>
    <expand macro="requirements" />
{%- if version_command %}
    <expand macro="version_command" />
{%- endif %}
{%- else %}
    <requirements>
{%- for requirement in requirements %}
        {{ requirement }}
{%- endfor %}
{%- for container in containers %}
        {{ container }}
{%- endfor %}
    </requirements>
{%- if version_command %}
    <version_command>{{ version_command }}</version_command>
{%- endif %}
{%- endif %}
    <command detect_errors="exit_code"><![CDATA[
{%- if command %}
        {{ command }}
{%- endif %}
{%- if not command and auto_commands %}
    TODO: executable name
{%- endif %}
{%- if auto_commands %}
        {{ auto_commands }}
{%- endif %}
{%- if not (command or auto_commands) %}
    TODO: Fill in command template.
{%- endif %}
    ]]></command>
    <inputs>
{%- for input in inputs %}
        {{ input }}
{%- endfor %}
{%- if auto_inputs %}
{{ auto_inputs }}
{%- endif %}
    </inputs>
    <outputs>
{%- for output in outputs %}
        {{ output }}
{%- endfor %}
    </outputs>
{%- if tests %}
    <tests>
{%- for test in tests %}
        <test>
{%- for param in test.params %}
            <param name="{{ param[0]}}" value="{{ param[1] }}"/>
{%- endfor %}
{%- for output in test.outputs %}
            <output name="{{ output[0] }}" file="{{ output[1] }}"/>
{%- endfor %}
        </test>
{%- endfor %}
    </tests>
{%- endif %}
    <help><![CDATA[
{%- if help %}
        {{ help }}
{%- endif %}
{%- if auto_help %}
        {{ auto_help }}
{%- endif %}
{%- if not (help or auto_help) %}
    TODO: Fill in help.
{%- endif %}
    ]]></help>
{%- if macros %}
    <expand macro="citations" />
{%- else %}
{%- if doi or bibtex_citations %}
    <citations>
{%- for single_doi in doi %}
        <citation type="doi">{{ single_doi }}</citation>
{%- endfor %}
{%- for bibtex_citation in bibtex_citations %}
        <citation type="bibtex">{{ bibtex_citation }}</citation>
{%- endfor %}
    </citations>
{%- endif %}
{%- endif %}
</tool>
"""

MACROS_TEMPLATE = """<macros>
    <xml name="requirements">
        <requirements>
{%- for requirement in requirements %}
        {{ requirement }}
{%- endfor %}
            <yield/>
{%- for container in containers %}
        {{ container }}
{%- endfor %}
        </requirements>
    </xml>
    <xml name="citations">
        <citations>
{%- for single_doi in doi %}
            <citation type="doi">{{ single_doi }}</citation>
{%- endfor %}
{%- for bibtex_citation in bibtex_citations %}
            <citation type="bibtex">{{ bibtex_citation }}</citation>
{%- endfor %}
            <yield />
        </citations>
    </xml>
{%- if version_command %}
    <xml name="version_command">
        <version_command>{{ version_command }}</version_command>
    </xml>
{%- endif %}
</macros>
"""

CWL_TEMPLATE = """#!/usr/bin/env cwl-runner
cwlVersion: '{{cwl_version}}'
class: CommandLineTool
id: "{{id}}"
label: "{{label}}"
{%- if containers or requirements %}
hints:
{%- for container in containers %}
  DockerRequirement:
    dockerPull: {{ container.image_id }}
{%- endfor %}
{%- if requirements %}
  SoftwareRequirement:
    packages:
{%- for requirement in requirements %}
    - package: {{ requirement.name }}
{%- if requirement.version %}
      version:
      - "{{ requirement.version }}"
{%- else %}
      version: []
{%- endif %}
{%- endfor %}
{%- endif %}
{%- endif %}
{%- if inputs or outputs %}
inputs:
{%- for input in inputs %}
  {{ input.id }}:
    type: {{ input.type }}
    doc: |
      TODO
    inputBinding:
      position: {{ input.position }}
{%- if input.prefix %}
      prefix: "{{input.prefix.prefix}}"
{%- if not input.prefix.separated %}
      separate: false
{%- endif %}
{%- endif %}
{%- endfor %}
{%- for output in outputs %}
{%- if output.require_filename %}
  {{ output.id }}:
    type: string
    doc: |
      Filename for output {{ output.id }}
    inputBinding:
      position: {{ output.position }}
{%- if output.prefix %}
      prefix: "{{output.prefix.prefix}}"
{%- if not output.prefix.separated %}
      separate: false
{%- endif %}
{%- endif %}
{%- endif %}
{%- endfor %}
{%- else %}
inputs: [] # TODO
{%- endif %}
{%- if outputs %}
outputs:
{%- for output in outputs %}
  {{ output.id }}:
    type: File
    outputBinding:
      glob: {{ output.glob }}
{%- endfor %}
{%- else %}
outputs: [] # TODO
{%- endif %}
{%- if base_command %}
baseCommand:
{%- for base_command_part in base_command %}
  - "{{ base_command_part}}"
{%- endfor %}
{%- else %}
baseCommand: []
{%- endif %}
{%- if arguments %}
arguments:
{%- for argument in arguments %}
  - valueFrom: "{{ argument.value }}"
    position: {{ argument.position }}
{%- if argument.prefix %}
      prefix: "{{argument.prefix.prefix}}"
{%- if not argument.prefix.separated %}
      separate: false
{%- endif %}
{%- endif %}
{%- endfor %}
{%- else %}
arguments: []
{%- endif %}
{%- if stdout %}
stdout: {{ stdout }}
{%- endif %}
doc: |
{%- if help %}
  {{ help|indent(2) }}
{%- else %}
   TODO: Fill in description.
{%- endif %}
"""

CWL_TEST_TEMPLATE = """
- doc: test generated from example command
  job: {{ job_filename }}
{%- if outputs %}
  outputs:
{%- for output in outputs %}
    {{ output.id }}:
      path: test-data/{{ output.example_value }}
{%- endfor %}
{%- else %}
  outputs: TODO
{%- endif %}
"""

CWL_JOB_TEMPLATE = """
{%- if inputs %}
{%- for input in inputs %}
{%- if input.type == "File" %}
{{ input.id }}:
  class: File
  path: test-data/{{ input.example_value }}
{%- else %}
  {{ input.id }}: {{ input.example_value }}
{%- endif %}
{%- endfor %}
{%- else %}
# TODO: Specify job input.
{}
{%- endif %}
"""


def build(**kwds):
    """Build up a :func:`ToolDescription` from supplid arguments."""
    if kwds.get("cwl"):
        builder = _build_cwl
    else:
        builder = _build_galaxy
    return builder(**kwds)


def _build_cwl(**kwds):
    _handle_help(kwds)
    _handle_requirements(kwds)
    assert len(kwds["containers"]) <= 1, kwds
    command_io = CommandIO(**kwds)
    render_kwds = {
        "cwl_version": DEFAULT_CWL_VERSION,
        "help": kwds.get("help", ""),
        "containers": kwds.get("containers", []),
        "requirements": kwds.get("requirements", []),
        "id": kwds.get("id"),
        "label": kwds.get("name"),
    }
    render_kwds.update(command_io.cwl_properties())

    contents = _render(render_kwds, template_str=CWL_TEMPLATE)
    tool_files = []
    test_files = []
    if kwds["test_case"]:
        sep = "-" if "-" in kwds.get("id") else "_"
        tests_path = f"{kwds.get('id')}{sep}tests.yml"
        job_path = f"{kwds.get('id')}{sep}job.yml"
        render_kwds["job_filename"] = job_path
        test_contents = _render(render_kwds, template_str=CWL_TEST_TEMPLATE)
        job_contents = _render(render_kwds, template_str=CWL_JOB_TEMPLATE)
        tool_files.append(ToolFile(tests_path, test_contents, "test"))
        tool_files.append(ToolFile(job_path, job_contents, "job"))
        for cwl_input in render_kwds["inputs"] or []:
            if cwl_input.type == "File" and cwl_input.example_value:
                test_files.append(cwl_input.example_value)

        for cwl_output in render_kwds["outputs"] or []:
            if cwl_output.example_value:
                test_files.append(cwl_output.example_value)

    return ToolDescription(contents, tool_files=tool_files, test_files=test_files)


def _build_galaxy(**kwds):
    # Test case to build up from supplied inputs and outputs, ultimately
    # ignored unless kwds["test_case"] is truthy.

    _handle_help(kwds)

    # process raw cite urls
    cite_urls = kwds.get("cite_url", [])
    del kwds["cite_url"]
    citations = list(map(UrlCitation, cite_urls))
    kwds["bibtex_citations"] = citations

    # handle requirements and containers
    _handle_requirements(kwds)

    command_io = CommandIO(**kwds)
    kwds["inputs"] = command_io.inputs
    kwds["outputs"] = command_io.outputs
    kwds["command"] = command_io.cheetah_template
    kwds["auto_inputs"] = command_io.auto_inputs
    kwds["auto_commands"] = command_io.auto_commands
    kwds["version_command"] = command_io.version_command
    test_case = command_io.test_case()

    # finally wrap up tests
    tests, test_files = _handle_tests(kwds, test_case)
    kwds["tests"] = tests

    # Render tool content from template.
    contents = _render(kwds)

    tool_files = []
    append_macro_file(tool_files, kwds)

    return ToolDescription(contents, tool_files=tool_files, test_files=test_files)


def append_macro_file(tool_files, kwds):
    macro_contents = None
    if kwds["macros"]:
        macro_contents = _render(kwds, MACROS_TEMPLATE)

        macros_file = "macros.xml"
        if not os.path.exists(macros_file):
            tool_files.append(ToolFile(macros_file, macro_contents, "macros"))

        io.info(REUSING_MACROS_MESSAGE)


class CommandIO:
    def __init__(self, **kwds):
        command = _find_command(kwds)
        cheetah_template = command

        # process raw inputs
        inputs = kwds.pop("input", [])
        inputs = list(map(Input, inputs or []))
        # alternatively process example inputs
        example_inputs = kwds.pop("example_input", [])
        for i, input_file in enumerate(example_inputs or []):
            name = f"input{i + 1}"
            inputs.append(Input(input_file, name=name, example=True))
            cheetah_template = _replace_file_in_command(cheetah_template, input_file, name)

        auto_inputs = None
        auto_commands = None
        auto_help = None
        version_command = None
        parser_path = kwds.get("autopygen", None)
        if parser_path is not None:
            parser = obtain_and_convert_parser(parser_path)

            if parser is not None:
                data_inputs = dict()
                reserved_names = set()
                name_map = dict()
                section_map = dict()

                generated_inputs, _, version_command_param = xml_from_decoy(
                    parser, data_inputs, reserved_names, name_map, section_map
                )

                auto_inputs = xml_to_string(generated_inputs, 8)
                # TODO make them useful  auto_outputs = xml_to_string(generated_outputs, 8)

                auto_commands = command_from_decoy(
                    parser, data_inputs, reserved_names, name_map, section_map, skip_default_namespace=True
                )

                if version_command_param:
                    version_command = f"[TODO exec name] {version_command_param.argument}"

                auto_help = parser.format_help()

        # handle raw outputs (from_work_dir ones) as well as named_outputs
        outputs = kwds.pop("output", [])
        outputs = list(map(Output, outputs or []))

        named_outputs = kwds.pop("named_output", [])
        for named_output in named_outputs or []:
            outputs.append(Output(name=named_output, example=False))

        # handle example outputs
        example_outputs = kwds.pop("example_output", [])
        for i, output_file in enumerate(example_outputs or []):
            name = f"output{i + 1}"
            from_path = output_file
            use_from_path = True
            if output_file in cheetah_template:
                # Actually found the file in the command, assume it can
                # be specified directly and skip from_work_dir.
                use_from_path = False
            output = Output(name=name, from_path=from_path, use_from_path=use_from_path, example=True)
            outputs.append(output)
            cheetah_template = _replace_file_in_command(cheetah_template, output_file, output.name)

        self.inputs = inputs
        self.outputs = outputs
        self.command = command
        self.auto_inputs = auto_inputs
        self.auto_commands = auto_commands
        self.auto_help = auto_help
        self.version_command = version_command
        self.cheetah_template = cheetah_template

    def example_input_names(self):
        for input in self.inputs:
            if input.example:
                yield input.input_description

    def example_output_names(self):
        for output in self.outputs:
            if output.example:
                yield output.example_path

    def cwl_lex_list(self):
        if not self.command:
            return []

        command_parts = shlex.split(self.command)
        parse_list = []

        input_count = 0
        output_count = 0

        index = 0

        prefixed_parts = []
        while index < len(command_parts):
            value = command_parts[index]
            eq_split = value.split("=")

            prefix = None
            if not _looks_like_start_of_prefix(index, command_parts):
                index += 1
            elif len(eq_split) == 2:
                prefix = Prefix(eq_split[0] + "=", False)
                value = eq_split[1]
                index += 1
            else:
                prefix = Prefix(value, True)
                value = command_parts[index + 1]
                index += 2
            prefixed_parts.append((prefix, value))

        for position, (prefix, value) in enumerate(prefixed_parts):
            if value in self.example_input_names():
                input_count += 1
                input = _CwlInput(
                    f"input{input_count}",
                    position,
                    prefix,
                    value,
                )
                parse_list.append(input)
            elif value in self.example_output_names():
                output_count += 1
                output = _CwlOutput(
                    f"output{output_count}",
                    position,
                    prefix,
                    value,
                )
                parse_list.append(output)
            elif prefix:
                param_id = prefix.prefix.lower().rstrip("=")
                type_ = param_type(value)
                input = _CwlInput(
                    param_id,
                    position,
                    prefix,
                    value,
                    type_=type_,
                )
                parse_list.append(input)
            else:
                part = _CwlCommandPart(value, position, prefix)
                parse_list.append(part)
        return parse_list

    def cwl_properties(self):
        base_command = []
        arguments = []
        inputs = []
        outputs = []

        lex_list = self.cwl_lex_list()

        index = 0
        while index < len(lex_list):
            token = lex_list[index]
            if isinstance(token, _CwlCommandPart):
                base_command.append(token.value)
            else:
                break
            index += 1

        while index < len(lex_list):
            token = lex_list[index]
            if token.is_token(">"):
                break
            token.position = index - len(base_command) + 1
            if isinstance(token, _CwlCommandPart):
                arguments.append(token)
            elif isinstance(token, _CwlInput):
                inputs.append(token)
            elif isinstance(token, _CwlOutput):
                token.glob = f"$(inputs.{token.id})"
                outputs.append(token)

            index += 1

        stdout = None
        if index < len(lex_list):
            token = lex_list[index]
            if token.is_token(">") and (index + 1) < len(lex_list):
                output_token = lex_list[index + 1]
                if not isinstance(output_token, _CwlOutput):
                    output_token = _CwlOutput("std_out", None)

                output_token.glob = "out"
                output_token.require_filename = False
                outputs.append(output_token)
                stdout = "out"
                index += 2
            else:
                io.warn("Example command too complex, you will need to build it up manually.")

        return {
            "inputs": inputs,
            "outputs": outputs,
            "arguments": arguments,
            "base_command": base_command,
            "stdout": stdout,
        }

    def test_case(self):
        test_case = TestCase()
        for input in self.inputs:
            if input.example:
                test_case.params.append((input.name, input.input_description))

        for output in self.outputs:
            if output.example:
                test_case.outputs.append((output.name, output.example_path))

        return test_case


def _looks_like_start_of_prefix(index, parts):
    value = parts[index]
    if len(value.split("=")) == 2:
        return True
    if index + 1 == len(parts):
        return False
    next_value = parts[index + 1]
    next_value_is_not_start = (len(value.split("=")) != 2) and next_value[0] not in ["-", ">", "<", "|"]
    return value.startswith("-") and next_value_is_not_start


Prefix = namedtuple("Prefix", ["prefix", "separated"])


class _CwlCommandPart:
    def __init__(self, value, position, prefix):
        self.value = value
        self.position = position
        self.prefix = prefix

    def is_token(self, value):
        return self.value == value


class _CwlInput:
    def __init__(self, id, position, prefix, example_value, type_="File"):
        self.id = id
        self.position = position
        self.prefix = prefix
        self.example_value = example_value
        self.type = type_

    def is_token(self, value):
        return False


class _CwlOutput:
    def __init__(self, id, position, prefix, example_value):
        self.id = id
        self.position = position
        self.prefix = prefix
        self.glob = None
        self.example_value = example_value
        self.require_filename = True

    def is_token(self, value):
        return False


def _render(kwds, template_str=TOOL_TEMPLATE):
    """Render template variables to generate the final tool."""
    return templates.render(template_str, **kwds)


def _replace_file_in_command(command, specified_file, name):
    """Replace example file with cheetah variable name.

    Be sure to single quote the name.
    """
    # TODO: check if the supplied variant was single quoted already.
    if f'"{specified_file}"' in command:
        # Sample command already wrapped filename in double quotes
        command = command.replace(f'"{specified_file}"', f"'${name}'")
    elif (f" {specified_file} ") in (" " + command + " "):
        # In case of spaces, best to wrap filename in double quotes
        command = command.replace(specified_file, f"'${name}'")
    else:
        command = command.replace(specified_file, f"${name}")
    return command


def _handle_help(kwds):
    """Convert supplied help parameters into a help variable for template.

    If help_text is supplied, use as is. If help is specified from a command,
    run the command and use that help text.
    """
    help_text = kwds.get("help_text")
    if not help_text:
        help_from_command = kwds.get("help_from_command")
        if help_from_command:
            p = subprocess.Popen(
                help_from_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True
            )
            help_text = p.communicate()[0]

    del kwds["help_text"]
    del kwds["help_from_command"]

    kwds["help"] = help_text


def _handle_tests(kwds, test_case):
    """Build tool test abstractions.

    Given state built up from handling rest of arguments (test_case) and
    supplied kwds - build tests for template and corresponding test files.
    """
    test_files = []
    if kwds["test_case"]:
        tests = [test_case]
        test_files.extend(map(lambda x: x[1], test_case.params))
        test_files.extend(map(lambda x: x[1], test_case.outputs))
    else:
        tests = []
    return tests, test_files


def _handle_requirements(kwds):
    """Build tool requirement abstractions.

    Convert requirements and containers specified from the command-line
    into abstract format for consumption by the template.
    """
    requirements = kwds["requirement"]
    del kwds["requirement"]
    requirements = list(map(Requirement, requirements or []))

    container = kwds["container"]
    del kwds["container"]
    containers = list(map(Container, container or []))

    kwds["requirements"] = requirements
    kwds["containers"] = containers


def _find_command(kwds):
    """Find base command from supplied arguments or just return `None`.

    If no such command was supplied (template will just replace this
    with a TODO item).
    """
    command = kwds.get("command")
    if not command:
        command = kwds.get("example_command", None)
        if command:
            del kwds["example_command"]
    return command


class UrlCitation:
    """Describe citation for tool."""

    def __init__(self, url):
        self.url = url

    def __str__(self):
        if "github.com" in self.url:
            return self._github_str()
        else:
            return self._url_str()

    def _github_str(self):
        url = self.url
        title = url.split("/")[-1]
        return f"""
@misc{{github{title},
  author = {{LastTODO, FirstTODO}},
  year = {{TODO}},
  title = {{{title}}},
  publisher = {{GitHub}},
  journal = {{GitHub repository}},
  url = {{{url}}},
}}"""

    def _url_str(self):
        url = self.url
        return f"""
@misc{{renameTODO,
  author = {{LastTODO, FirstTODO}},
  year = {{TODO}},
  title = {{TODO}},
  url = {{{url}}},
}}"""


class ToolDescription:
    """An description of the tool and related files to create."""

    def __init__(self, contents, tool_files=None, test_files=None):
        self.contents = contents
        self.tool_files = tool_files or []
        self.test_files = test_files or []


class ToolFile:
    def __init__(self, filename, contents, description):
        self.filename = filename
        self.contents = contents
        self.description = description


class Input:
    def __init__(self, input_description, name=None, example=False):
        parts = input_description.split(".")
        name = name or parts[0]
        if len(parts) > 0:
            datatype = ".".join(parts[1:])
        else:
            datatype = "data"

        self.input_description = input_description
        self.example = example
        self.name = name
        self.datatype = datatype

    def __str__(self):
        self.datatype = self.datatype.split(".")[-1]
        return f'<param type="data" name="{self.name}" format="{self.datatype}" />'


class Output:
    def __init__(self, from_path=None, name=None, use_from_path=False, example=False):
        if from_path:
            parts = from_path.split(".")
            name = name or parts[0]
            if len(parts) > 1:
                datatype = ".".join(parts[1:])
            else:
                datatype = "data"
        else:
            name = name
            datatype = "data"

        self.name = name
        self.datatype = datatype
        if use_from_path:
            self.from_path = from_path
        else:
            self.from_path = None
        self.example = example
        if example:
            self.example_path = from_path

    def __str__(self):
        if self.from_path:
            return self._from_path_str()
        else:
            return self._named_str()

    def _from_path_str(self):
        return f'<data name="{self.name}" format="{self.datatype}" from_work_dir="{self.from_path}" />'

    def _named_str(self):
        return f'<data name="{self.name}" format="{self.datatype}" />'


class Requirement:
    def __init__(self, requirement):
        parts = requirement.split("@", 1)
        if len(parts) > 1:
            name = parts[0]
            version = "@".join(parts[1:])
        else:
            name = parts[0]
            version = None
        self.name = name
        self.version = version

    def __str__(self):
        if self.version is not None:
            attrs = f' version="{self.version}"'
        else:
            attrs = ""
        return f'<requirement type="package"{attrs}>{self.name}</requirement>'


def param_type(value):
    if re.match(r"^\d+$", value):
        return "int"
    elif re.match(r"^\d+?\.\d+?$", value):
        return "float"
    else:
        return "string"


class Container:
    def __init__(self, image_id):
        self.type = "docker"
        self.image_id = image_id

    def __str__(self):
        return f'<container type="{self.type}">{self.image_id}</container>'


class TestCase:
    def __init__(self):
        self.params = []
        self.outputs = []


def write_tool_description(ctx, tool_description, **kwds):
    """Write a tool description to the file system guided by supplied CLI kwds."""
    tool_id = kwds.get("id")
    output = kwds.get("tool")
    if not output:
        extension = "cwl" if kwds.get("cwl") else "xml"
        output = f"{tool_id}.{extension}"
    if not io.can_write_to_path(output, **kwds):
        ctx.exit(1)

    io.write_file(output, tool_description.contents)
    io.info(f"Tool written to {output}")
    for tool_file in tool_description.tool_files:
        if tool_file.contents is None:
            continue

        path = tool_file.filename
        if not io.can_write_to_path(path, **kwds):
            ctx.exit(1)
        io.write_file(path, tool_file.contents)
        io.info(f"Tool {tool_file.description} written to {path}")

    macros = kwds["macros"]
    macros_file = "macros.xml"
    if macros and not os.path.exists(macros_file):
        io.write_file(macros_file, tool_description.macro_contents)
    elif macros:
        io.info(REUSING_MACROS_MESSAGE)
    if tool_description.test_files:
        if not os.path.exists("test-data"):
            io.info("No test-data directory, creating one.")
            os.makedirs("test-data")
        for test_file in tool_description.test_files:
            io.info(f"Copying test-file {test_file}")
            try:
                shutil.copy(test_file, "test-data")
            except Exception as e:
                io.info(f"Copy of {test_file} failed: {e}")


__all__ = (
    "build",
    "ToolDescription",
    "write_tool_description",
)
