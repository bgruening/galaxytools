# SPDX-License-Identifier: Apache-2.0
"""CWL Expression parsing."""
import asyncio
import copy
import inspect
import json
from collections.abc import Awaitable, MutableMapping
from typing import Any, Optional, Union, cast

from schema_salad.utils import json_dumps

from cwl_utils.errors import JavascriptException, SubstitutionError, WorkflowException
from cwl_utils.loghandler import _logger
from cwl_utils.sandboxjs import JSEngine, default_timeout, get_js_engine, param_re
from cwl_utils.types import CWLObjectType, CWLOutputType
from cwl_utils.utils import bytes2str_in_dicts


def _convert_dumper(string: str) -> str:
    return f"{json.dumps(string)} + "


def scanner(scan: str) -> Optional[tuple[int, int]]:
    """Find JS relevant punctuation in a string."""
    DEFAULT = 0
    DOLLAR = 1
    PAREN = 2
    BRACE = 3
    SINGLE_QUOTE = 4
    DOUBLE_QUOTE = 5
    BACKSLASH = 6

    i = 0
    stack = [DEFAULT]
    start = 0
    while i < len(scan):
        state = stack[-1]
        c = scan[i]

        if state == DEFAULT:
            if c == "$":
                stack.append(DOLLAR)
            elif c == "\\":
                stack.append(BACKSLASH)
        elif state == BACKSLASH:
            stack.pop()
            if stack[-1] == DEFAULT:
                return i - 1, i + 1
        elif state == DOLLAR:
            if c == "(":
                start = i - 1
                stack.append(PAREN)
            elif c == "{":
                start = i - 1
                stack.append(BRACE)
            else:
                stack.pop()
                i -= 1
        elif state == PAREN:
            if c == "(":
                stack.append(PAREN)
            elif c == ")":
                stack.pop()
                if stack[-1] == DOLLAR:
                    return start, i + 1
            elif c == "'":
                stack.append(SINGLE_QUOTE)
            elif c == '"':
                stack.append(DOUBLE_QUOTE)
        elif state == BRACE:
            if c == "{":
                stack.append(BRACE)
            elif c == "}":
                stack.pop()
                if stack[-1] == DOLLAR:
                    return start, i + 1
            elif c == "'":
                stack.append(SINGLE_QUOTE)
            elif c == '"':
                stack.append(DOUBLE_QUOTE)
        elif state == SINGLE_QUOTE:
            if c == "'":
                stack.pop()
            elif c == "\\":
                stack.append(BACKSLASH)
        elif state == DOUBLE_QUOTE:
            if c == '"':
                stack.pop()
            elif c == "\\":
                stack.append(BACKSLASH)
        i += 1

    if len(stack) > 1 and not (len(stack) == 2 and stack[1] in (BACKSLASH, DOLLAR)):
        raise SubstitutionError(
            "Substitution error, unfinished block starting at position {}: '{}' stack was {}".format(
                start, scan[start:], stack
            )
        )
    return None


def evaluator(
    js_engine: JSEngine,
    ex: str,
    obj: CWLObjectType,
    jslib: str,
    fullJS: bool,
    **kwargs: Any,
) -> Optional[CWLOutputType]:
    js_engine = js_engine or get_js_engine()
    expression_parse_exception = None

    if (match := param_re.match(ex)) is not None:
        first_symbol = match.group(1)
        first_symbol_end = match.end(1)

        if first_symbol_end + 1 == len(ex) and first_symbol == "null":
            return None
        try:
            if first_symbol not in obj:
                raise WorkflowException("%s is not defined" % first_symbol)

            if inspect.iscoroutinefunction(js_engine.regex_eval):
                return asyncio.get_event_loop().run_until_complete(
                    cast(
                        Awaitable[CWLOutputType],
                        js_engine.regex_eval(
                            first_symbol,
                            ex[first_symbol_end:-1],
                            cast(CWLOutputType, obj[first_symbol]),
                            **kwargs,
                        ),
                    )
                )
            else:
                return cast(
                    CWLOutputType,
                    js_engine.regex_eval(
                        first_symbol,
                        ex[first_symbol_end:-1],
                        cast(CWLOutputType, obj[first_symbol]),
                        **kwargs,
                    ),
                )
        except WorkflowException as werr:
            expression_parse_exception = werr
    if fullJS:
        if inspect.iscoroutinefunction(js_engine.eval):
            return asyncio.get_event_loop().run_until_complete(
                cast(Awaitable[CWLOutputType], js_engine.eval(ex, jslib, **kwargs))
            )
        else:
            return cast(CWLOutputType, js_engine.eval(ex, jslib, **kwargs))
    else:
        if expression_parse_exception is not None:
            raise JavascriptException(
                "Syntax error in parameter reference '%s': %s. This could be "
                "due to using Javascript code without specifying "
                "InlineJavascriptRequirement." % (ex[1:-1], expression_parse_exception)
            )
        else:
            raise JavascriptException(
                "Syntax error in parameter reference '%s'. This could be due "
                "to using Javascript code without specifying "
                "InlineJavascriptRequirement." % ex
            )


def interpolate(
    scan: str,
    rootvars: CWLObjectType,
    jslib: str = "",
    fullJS: bool = False,
    strip_whitespace: bool = True,
    escaping_behavior: int = 2,
    convert_to_expression: bool = False,
    js_engine: Optional[JSEngine] = None,
    **kwargs: Any,
) -> Optional[CWLOutputType]:
    """
    Interpolate and evaluate.

    Note: only call with convert_to_expression=True on CWL Expressions in $()
    form that need interpolation.
    """
    if strip_whitespace:
        scan = scan.strip()
    parts = []
    if convert_to_expression:
        dump = _convert_dumper
        parts.append("${return ")
    else:

        def dump(string: str) -> str:
            return string

    w = scanner(scan)
    while w:
        if convert_to_expression:
            parts.append(f'"{scan[0: w[0]]}" + ')  # noqa: B907
        else:
            parts.append(scan[0 : w[0]])

        if scan[w[0]] == "$":
            if not convert_to_expression:
                js_engine = js_engine or get_js_engine()
                e = evaluator(
                    js_engine, scan[w[0] + 1 : w[1]], rootvars, jslib, fullJS, **kwargs
                )
                if w[0] == 0 and w[1] == len(scan) and len(parts) <= 1:
                    return e

                leaf = json_dumps(e, sort_keys=True)
                if leaf[0] == '"':
                    leaf = json.loads(leaf)
                parts.append(leaf)
            else:
                parts.append(
                    "function(){var item ="
                    + scan[w[0] : w[1]][2:-1]
                    + '; if (typeof(item) === "string"){ return item; } '
                    "else { return JSON.stringify(item); }}() + "
                )
        elif scan[w[0]] == "\\":
            if escaping_behavior == 1:
                # Old behavior.  Just skip the next character.
                e = scan[w[1] - 1]
                parts.append(dump(e))
            elif escaping_behavior == 2:
                # Backslash quoting requires a three character lookahead.
                e = scan[w[0] : w[1] + 1]
                if e in ("\\$(", "\\${"):
                    # Suppress start of a parameter reference, drop the
                    # backslash.
                    parts.append(dump(e[1:]))
                    w = (w[0], w[1] + 1)
                elif e[1] == "\\":
                    # Double backslash, becomes a single backslash
                    parts.append(dump("\\"))
                else:
                    # Some other text, add it as-is (including the
                    # backslash) and resume scanning.
                    parts.append(dump(e[:2]))
            else:
                raise Exception("Unknown escaping behavior %s" % escaping_behavior)
        scan = scan[w[1] :]
        w = scanner(scan)
    if convert_to_expression:
        parts.append(f'"{scan}"')  # noqa: B907
        parts.append(";}")
    else:
        parts.append(scan)
    return "".join(parts)


def jshead(engine_config: list[str], rootvars: CWLObjectType) -> str:
    """Make sure all the byte strings are converted to str in `rootvars` dict."""
    return "\n".join(
        engine_config
        + [f"var {k} = {json_dumps(v, indent=4)};" for k, v in rootvars.items()]
    )


def needs_parsing(snippet: Any) -> bool:
    return isinstance(snippet, str) and ("$(" in snippet or "${" in snippet)


def do_eval(
    ex: Optional[CWLOutputType],
    jobinput: CWLObjectType,
    requirements: list[CWLObjectType],
    outdir: Optional[str],
    tmpdir: Optional[str],
    resources: dict[str, Union[float, int]],
    context: Optional[CWLOutputType] = None,
    timeout: float = default_timeout,
    strip_whitespace: bool = True,
    cwlVersion: str = "",
    **kwargs: Any,
) -> Optional[CWLOutputType]:
    """
    Evaluate the given CWL expression, in context.

    :param timeout: The maximum number of seconds to wait while executing.
    """
    runtime = cast(MutableMapping[str, Union[int, str, None]], copy.deepcopy(resources))
    runtime["tmpdir"] = tmpdir if tmpdir else None
    runtime["outdir"] = outdir if outdir else None

    rootvars = cast(
        CWLObjectType,
        bytes2str_in_dicts({"inputs": jobinput, "self": context, "runtime": runtime}),
    )

    if isinstance(ex, str) and needs_parsing(ex):
        fullJS = False
        jslib = ""
        for r in reversed(requirements):
            if r["class"] == "InlineJavascriptRequirement":
                fullJS = True
                jslib = jshead(cast(list[str], r.get("expressionLib", [])), rootvars)
                break

        try:
            return interpolate(
                ex,
                rootvars,
                timeout=timeout,
                fullJS=fullJS,
                jslib=jslib,
                strip_whitespace=strip_whitespace,
                escaping_behavior=(
                    1
                    if cwlVersion
                    in (
                        "v1.0",
                        "v1.1.0-dev1",
                        "v1.1",
                        "v1.2.0-dev1",
                        "v1.2.0-dev2",
                        "v1.2.0-dev3",
                    )
                    else 2
                ),
                **kwargs,
            )

        except Exception as e:
            _logger.exception(e)
            raise WorkflowException("Expression evaluation error:\n%s" % str(e)) from e
    else:
        return ex
