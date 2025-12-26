"""
Functions used to extract initialized ArgumentParser from target source code,
and to transform it into easily usable ParamInfo class containing
all the necessary info for 'inputs' element initialization
"""

import ast
import logging
import math
import re
import uuid
from argparse import ArgumentParser
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)

from lxml import etree

from planemo.autopygen.commands.command_utils import (
    create_element_with_body,
    transform_param_info,
)
from planemo.autopygen.param_info import (
    ParamDataType,
    ParamInfo,
    ParamTypeFlags,
)
from planemo.autopygen.source_file_parsing.decoy_parser import DecoyParser
from planemo.autopygen.source_file_parsing.local_module_parsing import handle_local_module_names
from planemo.autopygen.source_file_parsing.parser_discovery_and_init import get_parser_init_and_actions
from planemo.autopygen.source_file_parsing.parsing_commons import (
    create_module_tree_from_path,
    create_module_tree_from_str,
)
from planemo.autopygen.source_file_parsing.parsing_exceptions import ArgumentParsingDiscoveryError
from planemo.autopygen.source_file_parsing.unknown_names_discovery import initialize_variables_in_module
from planemo.autopygen.xml.xml_utils import (
    formatted_xml_elem,
    options,
    param,
    repeat,
)

DEFAULT_XML_INDENT = 4


def obtain_and_convert_parser(path: str) -> Optional[DecoyParser]:
    """
    Function parses python source code located at 'path',
    extracts initialized argument parser from the source file and returns it

    Parameters
    ----------
    path : str
      path to the source file containing argument parser init

    Returns
    -------
    Initialized argument parser, or None, in case error happened
    """
    try:
        tree = create_module_tree_from_path(path)
    except FileNotFoundError:
        logging.error("Input file not found")
        return None

    return obtain_parser(tree)


def obtain_and_convert_parser_from_str(text: str) -> Optional[ArgumentParser]:
    """
    Function parses python source code stored in text and
    extracts argument parser

    Parameters
    ----------
    text : str
      variable containing target source code

    Returns
    -------
    Initialized argument parser, or None, in case error happened
    """
    tree = create_module_tree_from_str(text)

    return obtain_parser(tree)


def obtain_parser(tree: ast.Module) -> Optional[DecoyParser]:
    try:
        actions, name, known_names = get_parser_init_and_actions(tree)

        actions, unknown_names = initialize_variables_in_module(tree, name, actions, known_names)

        result_module = handle_local_module_names(actions, unknown_names)
    except ArgumentParsingDiscoveryError as e:
        logging.error(e)
        return None

    ast.fix_missing_locations(result_module)
    compiled_module = compile(result_module, filename="<parser>", mode="exec")
    namespace: Dict[Any, Any] = {}
    try:
        exec(compiled_module, namespace)
    except Exception as e:
        print(e)
        logging.error("Parser couldn't be extracted")
        return None
    return namespace[name]


def _action_to_param(param_info: ParamInfo):
    if param_info.param_type.is_repeat:
        return repeat(param_info)
    elif param_info.param_type.is_selection:
        return options(param_info)

    return param(param_info)


def transform_section_name(section_name: str) -> str:
    return re.sub("[/\\-* ()]", "_", section_name).lower()


def generate_xml_from_section(
    section: DecoyParser.Section,
    data_inputs: Dict[str, str],
    reserved_names: Set[str],
    name_map: Dict[str, str],
    section_map: Dict[str, str],
    dont_wrap_in_section: bool = False,
) -> Tuple[List[etree._Element], List[Any], Optional[ParamInfo]]:
    sub_params = []
    sub_outputs = []
    version_command = None

    for action in section.actions:
        param_info = obtain_param_info(action, data_inputs, reserved_names, name_map, section_map)

        if param_info.param_type.is_help:
            continue

        if param_info.param_type.is_version:
            version_command = param_info
            continue

        sub_params.append(_action_to_param(param_info))

        if param_info.param_type.is_output:
            logging.warning("Outputs are not supported yet")
            continue

    sub_sections = []

    for subsection in (subsection for subsection in section.subsections if subsection.actions):
        inputs, outputs, sub_ver_comm = generate_xml_from_section(
            subsection, data_inputs, reserved_names, name_map, section_map
        )
        # we want to save the first version command found
        if sub_ver_comm is not None and version_command is None:
            version_command = sub_ver_comm

        sub_sections.extend(inputs)

        sub_outputs.extend(outputs)

    transformed_name = transform_section_name(section_map.get(section.name, section.name))

    result_inputs = [*sub_params, *sub_sections]

    if not dont_wrap_in_section:
        attrs = {
            "name": transformed_name,
            "title": section.name,
            "expanded": "true",
        }
        return [formatted_xml_elem("section", attrs, body=result_inputs)], sub_outputs, version_command

    return result_inputs, sub_outputs, version_command


def _command_recursion(
    section: DecoyParser.Section,
    data_inputs: Dict[str, str],
    reserved_names: Set[str],
    name_map: Dict[str, str],
    section_map: Dict[str, str],
    namespace: str,
    depth: int,
):
    for action in section.actions:
        param_info = obtain_param_info(action, data_inputs, reserved_names, name_map, section_map)

        if param_info.param_type.is_help or param_info.param_type.is_version:
            continue

        yield transform_param_info(param_info, namespace, depth)

    for subsection in section.subsections:
        sec_name = transform_section_name(section_map.get(subsection.name, subsection.name))
        separator = "." if namespace else ""
        child_namespace = f"{namespace}{separator}{sec_name}"

        yield from _command_recursion(
            subsection, data_inputs, reserved_names, name_map, section_map, child_namespace, depth + 1
        )


def _sub_parsers_conditionals(
    parsers: List[DecoyParser],
    index: int,
    data_inputs: Dict[str, str],
    reserved_names: Set[str],
    name_map: Dict[str, str],
    section_map: Dict[str, str],
    found_version_comm: bool = False,
) -> Tuple[List[Any], List[Any], Optional[ParamInfo]]:
    result = []
    result_outputs = []
    version_command = None
    conditional_options = []
    is_first = True
    for parser in parsers:
        conditional_options.append(
            formatted_xml_elem("option", {"value": parser.name, "selected": str(is_first).lower()}, text=parser.name)
        )
        if is_first:
            is_first = False

        inputs, outputs, version_comm = xml_from_decoy(parser, data_inputs, reserved_names, name_map, section_map)
        result_outputs.extend(outputs)
        result.append(formatted_xml_elem("when", {"value": parser.name}, body=inputs))

        if not found_version_comm or (version_command is None and version_comm is not None):
            version_command = version_comm

    conditional_param_definition = formatted_xml_elem(
        "param", {"name": "subparser_selector", "type": "select"}, body=conditional_options
    )
    result_inputs = formatted_xml_elem(
        "conditional", {"name": f"subparsers{index}"}, body=[conditional_param_definition, *result]
    )

    return [result_inputs], result_outputs, version_command


def xml_to_string(nodes: List[etree._Element], indent: int):
    def set_indent_and_convert(node: etree._Element) -> str:
        etree.indent(node, " " * DEFAULT_XML_INDENT)
        return etree.tostring(node, encoding="unicode", pretty_print=True)

    unified_str = "".join(map(set_indent_and_convert, nodes))
    return "\n".join(map(lambda line: " " * indent + line, unified_str.split("\n"))).rstrip(" ")


def xml_from_decoy(
    parser: DecoyParser,
    data_inputs: Dict[str, str],
    reserved_names: Set[str],
    name_map: Dict[str, str],
    section_map: Dict[str, str],
    dont_wrap_default_section: bool = True,
) -> Tuple[List[Any], List[Any], Optional[ParamInfo]]:
    input_params, outputs_list, version_command = generate_xml_from_section(
        parser.default_section, data_inputs, reserved_names, name_map, section_map, dont_wrap_default_section
    )
    inputs = [*input_params]

    outputs = [*outputs_list]
    for index, sub_parsers in enumerate(parser.sub_parsers):
        sub_parser_inputs, sub_parser_outputs, ver_comm = _sub_parsers_conditionals(
            sub_parsers.parsers, index, data_inputs, reserved_names, name_map, section_map, version_command is not None
        )

        inputs.extend(sub_parser_inputs)
        outputs.extend(sub_parser_outputs)

        if version_command is None and ver_comm is not None:
            version_command = ver_comm

    return inputs, outputs, version_command


def command_from_decoy(
    parser: DecoyParser,
    data_inputs: Dict[str, str],
    reserved_names: Set[str],
    name_map: Dict[str, str],
    section_map: Dict[str, str],
    depth: int = 0,
    starting_namespace: str = "",
    skip_default_namespace: bool = False,
) -> str:
    """
    The function generates commands from decoy parser. It requires name and section maps that have already been
    initialised and contain correct mapping
    """
    sec_name = starting_namespace
    if not skip_default_namespace:
        if starting_namespace:
            sec_name += "."

        sec_name += transform_section_name(section_map[parser.default_section.name])

    result = [
        "\n".join(
            _command_recursion(
                parser.default_section, data_inputs, reserved_names, name_map, section_map, sec_name, depth
            )
        )
    ]

    for index, sub_parsers in enumerate(parser.sub_parsers):
        subparsers_variable_name = f"subparsers{index}"
        subparsers_chosen_parser = f"${subparsers_variable_name}.subparser_selector"
        result.append(
            create_element_with_body(
                "if",
                subparsers_chosen_parser,
                [subparsers_chosen_parser],
                "Selected subparser",
                depth,
                body_indented=False,
            )
        )
        for parser in sub_parsers.parsers:
            conditional = create_element_with_body(
                "if",
                f'str({subparsers_chosen_parser}) == "{parser.name}"',
                body=[
                    command_from_decoy(
                        parser,
                        data_inputs,
                        reserved_names,
                        name_map,
                        section_map,
                        depth + 1,
                        subparsers_variable_name,
                        skip_default_namespace,
                    )
                ],
                comment=f"Conditional option {parser.name}",
                depth=0,
            )
            result.append(conditional)

    return "\n".join(result)


def obtain_param_info(
    action: DecoyParser.Action,
    data_inputs: Dict[str, str],
    reserved_names: Set[str],
    name_map: Dict[str, str],
    section_map: Dict[str, str],
) -> ParamInfo:
    name = action.argument.lstrip("-").replace("-", "_")
    argument = action.argument

    is_positional = name == argument

    nargs = _determine_nargs(action.kwargs.get("nargs", None))

    update_name_map(name, name_map, reserved_names)

    section = action.scope.name
    update_name_map(section, section_map, reserved_names)
    default_val = action.kwargs.get("default", None)

    help_ = action.kwargs.get("help", None)
    # value of required should be used, otherwise argparse assumes that all positional arguments are required
    optional = action.kwargs.get("required", not is_positional)

    param_type = _determine_param_type(action, nargs)
    type_ = _determine_type(action, data_inputs, name, action.kwargs.get("type", str), param_type)
    choices = action.kwargs.get("choices", None)
    custom_attributes = _determine_custom_attributes(param_type, type_, nargs)

    return ParamInfo(
        is_positional,
        param_type,
        type_,
        name_map[name],
        argument,
        name,
        section_map[section],
        section,
        default_val,
        custom_attributes,
        nargs,
        help_,
        optional,
        choices,
    )


def update_name_map(name: str, name_map: Dict[str, str], reserved_names: Set[str]):
    """
    Updates names that are equal to one of the values in 'reserved_names'
    Doesn't change names already present in the mapping

    Parameters
    ----------
    name :
    name_map :
    reserved_names :
    """
    old_name = name

    if old_name in name_map:
        return

    while name.lower() in reserved_names:
        name = name + str(uuid.uuid4())[:4]

    name_map[old_name] = name


def _determine_param_type(action: DecoyParser.Action, nargs: Union[int, float]):
    param_flags = ParamTypeFlags()
    # TODO add output
    param_flags.is_selection = action.kwargs.get("choices", None) is not None
    param_flags.is_repeat = action.action in ["APPEND", "APPEND_CONST", "COUNT"] or nargs > 1
    param_flags.is_flag = action.action in ["STORE_TRUE", "STORE_CONST", "STORE_FALSE", "APPEND_CONST", "COUNT"]
    param_flags.is_extend = action.action == "EXTEND"
    param_flags.is_help = action.action == "HELP"
    param_flags.is_version = action.action == "VERSION"

    return param_flags


def _determine_type(
    action: DecoyParser.Action, data_inputs: Dict[str, str], name: str, type_, param_type_flags: ParamTypeFlags
) -> ParamDataType:
    if action.kwargs.get("choices", []):
        return ParamDataType.SELECT

    if type_ is None or type_ is bool or param_type_flags.is_flag:
        return ParamDataType.BOOLEAN
    elif type_ is int:
        return ParamDataType.INTEGER
    elif type_ is float:
        return ParamDataType.FLOAT
    elif type_ is str:
        if name in data_inputs:
            return ParamDataType.DATA
        else:
            return ParamDataType.TEXT

    return ParamDataType.UNDEFINED


def _determine_nargs(nargs: Union[str, int, None]) -> Union[float, int]:
    if isinstance(nargs, str):
        if nargs == "?":
            return 1
        return math.inf
    if nargs is None:
        return 0
    return int(nargs)


def _determine_custom_attributes(
    param_type: ParamTypeFlags, type_: ParamDataType, nargs: Union[float, int]
) -> Dict[str, str]:
    flags: Dict[str, str] = dict()

    if type_ == ParamDataType.DATA and nargs > 1:
        flags["multiple"] = "true"

    if param_type.is_repeat and 1 < nargs < math.inf:
        flags["max"] = str(nargs)

    return flags
