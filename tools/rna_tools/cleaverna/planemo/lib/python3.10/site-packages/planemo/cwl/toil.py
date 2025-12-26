import json
import tempfile

from planemo.deps import ensure_dependency_resolvers_conf_configured
from planemo.io import (
    error,
    real_io,
)
from planemo.runnable import (
    ErrorRunResponse,
    Runnable,
)
from .run import (
    CwlToolRunResponse,
    JSON_PARSE_ERROR_MESSAGE,
)

TOIL_REQUIRED_MESSAGE = "This functionality requires Toil, please install with 'pip install toil'"


def run_toil(ctx, runnable: "Runnable", job_path: str, **kwds):
    """Translate planemo kwds to cwltool kwds and run cwltool main function."""
    cwltoil = _ensure_toil_available()

    args = []
    if not ctx.verbose:
        args.append("--quiet")
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
    ctx.vlog(f"Calling cwltoil with arguments {args}")
    with tempfile.NamedTemporaryFile("w") as tmp_stdout:
        # cwltool passes sys.stderr to subprocess.Popen - ensure it has
        # and actual fileno.
        with real_io():
            ret_code = cwltoil.main(args, stdout=tmp_stdout)
        tmp_stdout.flush()
        with open(tmp_stdout.name) as stdout_f:
            try:
                result = json.load(stdout_f)
            except ValueError:
                message = JSON_PARSE_ERROR_MESSAGE % (
                    open(tmp_stdout.name).read(),
                    tmp_stdout.name,
                )
                error(message)
                raise Exception(message)

        if ret_code != 0:
            return ErrorRunResponse("Error running Toil")
        outputs = result
    return CwlToolRunResponse(
        runnable,
        "",
        outputs=outputs,
    )


def _ensure_toil_available():
    try:
        from toil.cwl import cwltoil

        return cwltoil
    except ImportError:
        raise Exception(TOIL_REQUIRED_MESSAGE)
