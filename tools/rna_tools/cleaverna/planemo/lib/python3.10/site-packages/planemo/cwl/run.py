"""Module defines a planemo abstraction around running cwltool.

cwltool is an executable Python script and library mostly maintained by
Peter Amstutz and serves the reference implementation for the CWL.
It can be found at https://github.com/common-workflow-language/cwltool,
"""

import json
import tempfile
from typing import (
    Any,
    Dict,
    Optional,
    TYPE_CHECKING,
    Union,
)

from galaxy.tool_util.cwl.cwltool_deps import (
    ensure_cwltool_available,
    main,
)

from planemo.deps import ensure_dependency_resolvers_conf_configured
from planemo.io import (
    error,
    real_io,
)
from planemo.runnable import (
    ErrorRunResponse,
    Runnable,
    SuccessfulRunResponse,
)

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

JSON_PARSE_ERROR_MESSAGE = "Failed to parse JSON from cwltool output [%s] in file [%s]. cwltool logs [%s]."


class CwlToolRunResponse(SuccessfulRunResponse):
    """Describe the resut of a cwltool invocation."""

    def __init__(
        self,
        runnable: "Runnable",
        log: str,
        outputs: Optional[Dict[str, Any]] = None,
    ) -> None:
        super().__init__(runnable=runnable)
        self._log = log
        self._outputs = outputs

    @property
    def log(self) -> str:
        return self._log

    @property
    def job_info(self) -> None:
        return None

    @property
    def invocation_details(self) -> None:
        return None

    @property
    def outputs_dict(self):
        return self._outputs


def run_cwltool(
    ctx: "PlanemoCliContext", runnable: "Runnable", job_path: str, **kwds
) -> Union[ErrorRunResponse, CwlToolRunResponse]:
    """Translate planemo kwds to cwltool kwds and run cwltool main function."""
    ensure_cwltool_available()

    args = []
    if ctx.verbose:
        args.append("--verbose")
    output_directory = kwds.get("output_directory", None)
    if output_directory:
        args.append("--outdir")
        args.append(output_directory)
    if kwds.get("no_container", False):
        args.append("--no-container")
        ensure_dependency_resolvers_conf_configured(ctx, kwds)
        args.append("--beta-dependency-resolvers-configuration")
        args.append(kwds["dependency_resolvers_config_file"])
    if kwds.get("mulled_containers"):
        args.append("--beta-use-biocontainers")

    if kwds.get("non_strict_cwl", False):
        args.append("--non-strict")

    args.extend([runnable.path, job_path])
    ctx.vlog("Calling cwltool with arguments %s" % args)
    with tempfile.NamedTemporaryFile("w") as tmp_stdout, tempfile.NamedTemporaryFile("w") as tmp_stderr:
        # cwltool passes sys.stderr to subprocess.Popen - ensure it has
        # and actual fileno.
        with real_io():
            ret_code = main.main(
                args,
                stdout=tmp_stdout,
                stderr=tmp_stderr,
            )
        tmp_stdout.flush()
        tmp_stderr.flush()
        with open(tmp_stderr.name) as stderr_f:
            log = stderr_f.read()
            ctx.vlog("cwltool log output [%s]" % log)
        with open(tmp_stdout.name) as stdout_f:
            try:
                result = json.load(stdout_f)
            except ValueError:
                message = JSON_PARSE_ERROR_MESSAGE % (open(tmp_stdout.name).read(), tmp_stdout.name, log)
                error(message)
                raise Exception(message)

        if ret_code != 0:
            return ErrorRunResponse("Error running cwltool", log=log)
        outputs = result
    return CwlToolRunResponse(
        runnable,
        log,
        outputs=outputs,
    )


__all__ = ("run_cwltool",)
