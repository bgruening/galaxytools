"""
Module used for creation and manipulation of template elements
are parts of template, separated from the rest by comments
with specific structure
Example of the comments:

    ## foo definition
    ... block itself ...
    ## end foo definition

Elements can be nested
"""

from typing import List

from planemo.autopygen.param_info import (
    ParamDataType,
    ParamInfo,
)

SPACE = " "
DEFAULT_INDENT = 4
ADD_COMMENTS_BY_DEFAULT = False


class DefinitionNotFoundException(Exception):
    """
    Exception raised if part of template definition cannot be found
    """

    pass


def create_flag(
    variable: str, comment: str, depth: int, indent=DEFAULT_INDENT, add_comment: bool = ADD_COMMENTS_BY_DEFAULT
) -> str:
    """
    Function used to create a flag definition, wrapped in a comment

    Parameters
    ----------
    variable : str
        name of variable, containing $ at the beginning
    comment : str
        wrapping comment
    depth : int
      integer, used to set the depth of the current element.
      This value is used to indent the block properly
    indent : int
      default value for size of the block indent
    add_comment : bool
      option that enables or disables formatting comments
    """
    result = f"{depth * indent * SPACE}{variable}\n"

    if not add_comment:
        return result

    return f"{depth * indent * SPACE}## FLAG {comment}\n" f"{result}" f"{depth * indent * SPACE}## end FLAG {comment}\n"


def create_element_with_body(
    kind: str,
    head: str,
    body: List[str],
    comment: str,
    depth: int,
    indent: int = DEFAULT_INDENT,
    body_indented: bool = True,
    add_comment: bool = ADD_COMMENTS_BY_DEFAULT,
) -> str:
    """
    Function used to create block of template, like if or loop

    Parameters
    ----------
    kind : str
      string defining what kind of element is created, for example if or for
      (loop)
    head : str
      body of block header, for example predicate of condition, or the
      body of loop
    body : str
      body of the block, can be another element
    comment : str
      comment, used to set the start and end of the block
    depth : int
      integer, used to set the depth of the current element.
      This value is used to indent the block properly
    indent : int
      default value for size of the block indent
    body_indented : bool
      option that define whether body is already correctly indented
    add_comment : bool
      option that enables or disables formatting comments

    Returns
    -------
    string containing the created template element
    """
    result = []
    if add_comment:
        result.append(f"{depth * indent * SPACE}## {comment}\n")

    result.append(f"{depth * indent * SPACE}#{kind} {head}:\n")

    body_indent = ""
    if body:
        if not body_indented:
            body_indent = (depth + 1) * indent * SPACE
            body[0] = f"{body_indent}{body[0]}"

        translated_body = ("\n" + body_indent).join(body)
        if translated_body[-1] != "\n":
            translated_body += "\n"

        result.append(translated_body)

    result.append(f"{depth * indent * SPACE}#end {kind}\n")

    if add_comment:
        result += f"{depth * indent * SPACE}## end {comment}\n"

    return "".join(result)


def transform_param_info(info: ParamInfo, namespace: str, depth: int):
    if info.param_type.is_help or info.param_type.is_version:
        raise ParamTypeNotSupported("Transformation for these param types are not supported")

    name = info.name
    separator = "." if namespace else ""
    variable = f"${namespace}{separator}{name}"
    if not info.param_type.is_repeat:
        if info.param_type.is_flag:
            return create_flag(variable, f"{name} definition", depth)
        else:
            body_expression = create_body_expression(info, variable, depth + 1)
            return create_element_with_body("if", variable, [body_expression], f"{name} definition", depth)

    iteration_var = "$item"
    if info.param_type.is_flag:
        param = create_flag(iteration_var, f"{name} definition", depth + 1)
    else:
        body_expression = create_body_expression(info, iteration_var, depth + 2)
        param = create_element_with_body("if", iteration_var, [body_expression], f"{name} definition", depth + 1)

    head_expression = f"{iteration_var} in ${namespace}{separator}{info.name}"

    return create_element_with_body("for", head_expression, [param], f"{info.name} definition", depth)


# TODO generating command like this assumes that parameters are added to argparse in the right order.
# Separating arguments into positional and non-positional is necessary,
# otherwise the command generator will not work correctly
def create_body_expression(info: ParamInfo, variable: str, depth: int, indentation: int = DEFAULT_INDENT) -> str:
    stripped_arg = info.argument.lstrip("-")
    str_indent = SPACE * depth * indentation

    wrapped_variable = variable
    if info.type == ParamDataType.DATA or info.type == ParamDataType.TEXT:
        wrapped_variable = f"'{wrapped_variable}'"

    if stripped_arg == info.argument:
        return f"{str_indent}{variable}"
    return f"{str_indent}{info.argument} {wrapped_variable}"


class ParamTypeNotSupported(Exception):
    pass
