"""Click definitions for various shared options and arguments."""

import functools
import os

import click
from galaxy.tool_util.deps import (
    docker_util,
    singularity_util,
)
from galaxy.tool_util.verify.interactor import DEFAULT_TOOL_TEST_WAIT
from gxjobconfinit.types import Runner

from .config import planemo_option


def force_option(what="files"):
    """Annotate click command as consume the -f/--force option."""
    return planemo_option(
        "-f",
        "--force",
        is_flag=True,
        help="Overwrite existing %s if present." % what,
    )


def open_file_option():
    return planemo_option(
        "--open",
        is_flag=True,
        help="Open the file in your default editor after creation.",
    )


def skip_venv_option():
    """Annotate click command as consume the --skip_venv option."""
    return planemo_option(
        "--skip_venv",
        is_flag=True,
        help=(
            "Do not create or source a virtualenv environment for Galaxy, "
            "this should be used to preserve an externally configured "
            "virtual environment or conda environment."
        ),
    )


def skip_client_build_option():
    """Annotate click command as consume the --skip_client_build option."""
    return planemo_option(
        "--skip_client_build",
        "galaxy_skip_client_build",
        is_flag=True,
        default=False,
        help=("Do not build Galaxy client when serving Galaxy."),
    )


def install_prebuilt_client_option():
    return planemo_option(
        "--install_prebuilt_client/--no_install_prebuilt_client",
        "galaxy_install_prebuilt_client",
        is_flag=True,
        default=True,
        help=("Install a pre-built client from npm. Turn this off you need access to visualizations."),
    )


def run_engine_option():
    """Annotate click command as consume the --engine option."""
    return planemo_option(
        "--engine",
        type=click.Choice(["galaxy", "docker_galaxy", "cwltool", "toil", "external_galaxy"]),
        default=None,
        use_global_config=True,
        help=(
            "Select an engine to run or test artifacts such as tools "
            "and workflows. Defaults to a local Galaxy, but running Galaxy within "
            "a Docker container or the CWL reference implementation 'cwltool' and "
            "'toil' be selected."
        ),
    )


def non_strict_cwl_option():
    """Annotate click command as consume the --non_strict_cwl option."""
    return planemo_option(
        "--non_strict_cwl",
        default=False,
        is_flag=True,
        help="Disable strict validation of CWL.",
    )


def serve_engine_option():
    """Annotate click command as consume the --engine option.

    This variant of the engine command is restricted to engines that can serve Galaxy
    servers.
    """
    return planemo_option(
        "--engine",
        type=click.Choice(["galaxy", "docker_galaxy", "external_galaxy"]),
        default="galaxy",
        use_global_config=True,
        use_env_var=True,
        help=(
            "Select an engine to serve artifacts such as tools "
            "and workflows. Defaults to a local Galaxy, but running Galaxy within "
            "a Docker container."
        ),
    )


def ignore_dependency_problems_option():
    """Annotate click command as consume the --ignore_dependency_problems option."""
    return planemo_option(
        "--ignore_dependency_problems",
        is_flag=True,
        default=False,
        use_global_config=True,
        help=(
            "When installing shed repositories for workflows, ignore dependency issues. "
            "These likely indicate a problem but in some cases may not prevent a workflow "
            "from successfully executing."
        ),
    )


def cwltool_no_container_option():
    """Annotate click command as consume the --no_container option.

    This option is for the CWL CLI runner interface.
    """
    return planemo_option(
        "--no-container",
        "--no_container",
        is_flag=True,
        default=False,
        use_global_config=True,
        help=("If cwltool engine is used, disable Docker container usage."),
    )


def test_data_option():
    """Annotate click command as consume the --test_data option."""
    return planemo_option(
        "--test_data",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help="test-data directory to for specified tool(s).",
    )


def extra_tools_option():
    return planemo_option(
        "--extra_tools",
        type=click.Path(exists=True, file_okay=True, dir_okay=True, resolve_path=True),
        multiple=True,
        help=(
            "Extra tool sources to include in Galaxy's tool panel (file or "
            "directory). These will not be linted/tested/etc... but they "
            "will be available to workflows and for interactive use."
        ),
    )


def tool_data_table_option():
    return planemo_option(
        "--tool_data_table",
        type=click.Path(exists=True, file_okay=True, resolve_path=True),
        help="tool_data_table_conf.xml file to for specified tool(s).",
        multiple=True,
    )


def galaxy_email_option():
    return planemo_option(
        "--galaxy_email",
        type=str,
        default="planemo@galaxyproject.org",
        use_global_config=True,
        use_env_var=True,
        help="E-mail address to use when launching single-user Galaxy server.",
    )


def galaxy_python_version():
    return planemo_option(
        "--galaxy_python_version",
        use_global_config=True,
        default=None,
        type=click.Choice(["3", "3.8", "3.9", "3.10", "3.11", "3.12"]),
        help="Python version to start Galaxy under",
    )


def galaxy_root_option():
    return planemo_option(
        "--galaxy_root",
        use_global_config=True,
        extra_global_config_vars=["galaxy_root"],
        use_env_var=True,
        type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
        help="Root of development galaxy directory to execute command with.",
    )


def galaxy_version_option():
    return planemo_option(
        "--galaxy_version",
        type=str,
        default="24.2",
        help="Version of Galaxy to target for configuration (default 24.2).",
    )


def tpv_option():
    return planemo_option(
        "--tpv/--no_tpv",
        is_flag=True,
        default=False,
        help=(
            "Include TPV (Total Perspective Vortex) configuration and shared usegalaxy* "
            "database of tool cores and memory for allocation purposes. "
        ),
    )


def galaxy_cwl_root_option():
    return planemo_option(
        "--cwl_galaxy_root",
        use_global_config=True,
        extra_global_config_vars=["cwl_galaxy_root"],
        use_env_var=True,
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help=(
            "Root of development galaxy directory to execute command with"
            " (must be branch of Galaxy with CWL support, this option"
            " is experimental and will be replaced with --galaxy_root when"
            " and if CWL support is merged into Galaxy."
        ),
    )


def galaxy_port_option():
    return planemo_option(
        "--port",
        type=int,
        default="9090",
        use_global_config=True,
        help="Port to serve Galaxy on (default is 9090).",
    )


def galaxy_host_option():
    return planemo_option(
        "--host",
        type=str,
        default="127.0.0.1",
        use_global_config=True,
        help=(
            "Host to bind Galaxy to. Default is 127.0.0.1 that is "
            "restricted to localhost connections for security reasons "
            "set to 0.0.0.0 to bind Galaxy to all ports including "
            "potentially publicly accessible ones."
        ),
    )


def dependency_resolvers_option():
    return planemo_option(
        "--dependency_resolvers_config_file",
        type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
        use_global_config=True,
        help="Dependency resolver configuration for Galaxy to target.",
    )


def enable_cwl_option():
    return planemo_option(
        "--cwl",
        is_flag=True,
        help=(
            "Configure Galaxy for use with CWL tool."
            " (this option is experimental and will be replaced when"
            " and if CWL support is merged into Galaxy)."
        ),
    )


def build_cwl_option():
    return planemo_option(
        "--cwl",
        is_flag=True,
        help="Build a CWL tool instead of a Galaxy tool.",
    )


def run_history_tags_option():
    return planemo_option(
        "tags",
        "--tags",
        type=str,
        default=None,
        help="Comma-separated list of tags to add to the created history.",
    )


def run_output_directory_option():
    return planemo_option(
        "output_directory",
        "--output_directory",
        "--outdir",
        type=click.Path(
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
        default=None,
        help=("Where to store outputs of a 'run' task."),
    )


def run_output_json_option():
    return planemo_option(
        "output_json",
        "--output_json",
        type=click.Path(
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
        default=None,
        help=("Where to store JSON dictionary describing outputs of a 'run' task."),
    )


def run_download_outputs_option():
    return planemo_option(
        "--download_outputs/--no_download_outputs",
        is_flag=True,
        default=False,
        help=(
            "After tool or workflow runs are complete, download "
            "the output files to the location specified by --output_directory. "
        ),
    )


def run_export_option():
    return planemo_option(
        "--export_invocation",
        help="Export workflow invocation as archive to specified path.",
        type=click.Path(),
        default=None,
    )


def publish_dockstore_option():
    return planemo_option(
        "--publish/--no_publish",
        is_flag=True,
        default=True,
        help="Set publish attribute to true in .dockstore.yml file",
    )


def no_dependency_resolution():
    return planemo_option(
        "--no_dependency_resolution",
        use_global_config=True,
        is_flag=True,
        help="Configure Galaxy with no dependency resolvers.",
    )


def brew_dependency_resolution():
    return planemo_option(
        "--brew_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use plain brew dependency resolution.",
    )


def conda_dependency_resolution():
    return planemo_option(
        "--conda_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use only conda for dependency resolution.",
    )


def shed_dependency_resolution():
    return planemo_option(
        "--shed_dependency_resolution",
        is_flag=True,
        help=("Configure Galaxy to use brewed Tool Shed dependency resolution."),
    )


def file_path_option():
    return planemo_option(
        "--file_path",
        type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
        help="Location for files created by Galaxy (e.g. database/files).",
        default=None,
        use_global_config=True,
    )


def database_connection_option():
    return planemo_option(
        "--database_connection",
        type=str,
        help="Database connection string to use for Galaxy.",
        default=None,
        use_global_config=True,
    )


def shed_tools_conf_option():
    return planemo_option(
        "--shed_tool_conf",
        type=str,
        help="Location of shed tools conf file for Galaxy.",
        default=None,
        use_global_config=True,
    )


def shed_tools_directory_option():
    return planemo_option(
        "--shed_tool_path",
        type=str,
        help="Location of shed tools directory for Galaxy.",
        default=None,
        use_global_config=True,
    )


def tool_dependency_dir_option():
    return planemo_option(
        "--tool_dependency_dir",
        type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
        default=None,
        use_global_config=True,
        help="Tool dependency dir for Galaxy to target.",
    )


def job_config_option():
    return planemo_option(
        "--job_config_file",
        type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
        help="Job configuration file for Galaxy to target.",
        default=None,
        use_global_config=True,
    )


class EnumType(click.Choice):
    def __init__(self, enum):
        self._enum = enum
        super().__init__([e.value for e in enum])

    def convert(self, value, param, ctx):
        return self._enum(super().convert(value, param, ctx))


def runner_target_option():
    return planemo_option(
        "--runner",
        type=EnumType(Runner),
        help="Galaxy runner (e.g. DRM) to target.",
        default=Runner.LOCAL,
        use_global_config=False,
    )


def tool_data_path_option():
    return planemo_option(
        "--tool_data_path",
        type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
        help="Directory where data used by tools is located. Required if tests are run in docker and should make use of external reference data.",
        default=None,
        use_global_config=True,
    )


def mulled_containers_option():
    return planemo_option(
        "mulled_containers",
        "--mulled_containers",
        "--biocontainers",
        is_flag=True,
        help="Test tools against mulled containers (forces --docker). Disables conda resolution unless any conda option has been set explicitly.",
    )


def galaxy_startup_timeout_option():
    return planemo_option(
        "--galaxy_startup_timeout",
        type=click.IntRange(1),
        default=900,
        help="Wait for galaxy to start before assuming Galaxy did not start.",
    )


def install_galaxy_option():
    return planemo_option(
        "--install_galaxy", is_flag=True, help="Download and configure a disposable copy of Galaxy from github."
    )


def docker_galaxy_image_option():
    return planemo_option(
        "--docker_galaxy_image",
        default="quay.io/bgruening/galaxy",
        use_global_config=True,
        help=(
            "Docker image identifier for docker-galaxy-flavor used if "
            "engine type is specified as ``docker-galaxy``. Defaults to "
            "quay.io/bgruening/galaxy."
        ),
    )


def docker_extra_volume_option():
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    )
    return planemo_option(
        "--docker_extra_volume",
        type=arg_type,
        default=None,
        use_global_config=True,
        multiple=True,
        help=("Extra path to mount if --engine docker or `--biocontainers` or `--docker`."),
    )


def singularity_extra_volume_option():
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    )
    return planemo_option(
        "--singularity_extra_volume",
        type=arg_type,
        default=None,
        use_global_config=True,
        multiple=True,
        help=("Extra path to mount if --engine docker or `--biocontainers` or `--singularity`."),
    )


def galaxy_url_option(required: bool = False):
    return planemo_option(
        "--galaxy_url",
        use_global_config=True,
        extra_global_config_vars=["galaxy_url"],
        use_env_var=True,
        required=required,
        type=str,
        help="Remote Galaxy URL to use with external Galaxy engine.",
    )


def galaxy_admin_key_option():
    return planemo_option(
        "--galaxy_admin_key",
        use_global_config=True,
        extra_global_config_vars=["admin_key"],
        use_env_var=True,
        type=str,
        help="Admin key to use with external Galaxy engine.",
    )


def galaxy_user_key_option(required: bool = False):
    return planemo_option(
        "--galaxy_user_key",
        use_global_config=True,
        extra_global_config_vars=["admin_key"],
        use_env_var=True,
        required=required,
        type=str,
        help="User key to use with external Galaxy engine.",
    )


def history_name():
    return planemo_option(
        "--history_name",
        type=str,
        help="Name to give a Galaxy history, if one is created.",
    )


def history_id():
    return planemo_option(
        "--history_id",
        type=str,
        help="Send the results of the run to the history with the provided ID. A history with this ID must exist.",
    )


def no_cache_galaxy_option():
    return planemo_option(
        "--no_cache_galaxy",
        is_flag=True,
        help=(
            "Skip caching of Galaxy source and dependencies obtained with "
            "--install_galaxy. Not caching this results in faster "
            "downloads (no git) - so is better on throw away instances such "
            "with TravisCI. "
        ),
    )


def galaxy_branch_option():
    return planemo_option(
        "--galaxy_branch",
        default=None,
        use_global_config=True,
        use_env_var=True,
        help=("Branch of Galaxy to target (defaults to master) if a Galaxy root isn't specified."),
    )


def galaxy_source_option():
    return planemo_option(
        "--galaxy_source",
        default=None,
        use_global_config=True,
        help=(
            "Git source of Galaxy to target (defaults to the official "
            "galaxyproject github source if a Galaxy root isn't "
            "specified."
        ),
    )


def skip_install_option():
    return planemo_option(
        "--skip_install", is_flag=True, help="Skip installation - only source requirements already available."
    )


def brew_option():
    return planemo_option(
        "--brew",
        use_global_config=True,
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Homebrew 'brew' executable to use.",
    )


def conda_prefix_option():
    return planemo_option(
        "--conda_prefix",
        use_global_config=True,
        use_env_var=True,
        type=click.Path(file_okay=False, dir_okay=True),
        help="Conda prefix to use for conda dependency commands.",
    )


def conda_exec_option():
    return planemo_option(
        "--conda_exec",
        use_global_config=True,
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Location of conda executable.",
    )


def conda_use_local_option():
    return planemo_option(
        "--conda_use_local", is_flag=True, help="Use locally built packages while building Conda environments."
    )


def conda_ensure_channels_option():
    return planemo_option(
        "conda_ensure_channels",
        "--conda_channels",
        "--conda_ensure_channels",
        type=str,
        use_global_config=True,
        use_env_var=True,
        help=("Ensure conda is configured with specified comma separated list of channels."),
        default="conda-forge,bioconda",
    )


def conda_auto_install_option():
    return planemo_option(
        "--conda_auto_install/--no_conda_auto_install",
        is_flag=True,
        default=True,
        help=("Conda dependency resolution for Galaxy will attempt to install requested but missing packages."),
    )


def conda_auto_init_option():
    return planemo_option(
        "--conda_auto_init/--no_conda_auto_init",
        is_flag=True,
        default=True,
        help=(
            "Conda dependency resolution for Galaxy will auto install "
            "conda itself using miniforge if not availabe on conda_prefix."
        ),
    )


def conda_global_option():
    return planemo_option(
        "--global",
        is_flag=True,
        default=False,
        help=(
            "Install Conda dependencies globally instead of in requirement specific "
            "environments packaged for tools. If the Conda bin directory is on your "
            "PATH, tools may still use binaries but this is more designed for "
            "interactive testing and debugging."
        ),
    )


def simultaneous_upload_option():
    return planemo_option(
        "--simultaneous_uploads/--no_simultaneous_uploads",
        is_flag=True,
        use_global_config=True,
        default=False,
        help=(
            "When uploading files to Galaxy for tool or workflow tests or runs, "
            "upload multiple files simultaneously without waiting for the previous "
            "file upload to complete."
        ),
    )


def check_uploads_ok_option():
    return planemo_option(
        "--check_uploads_ok/--no_check_uploads_ok",
        is_flag=True,
        default=True,
        help=(
            "When uploading files to Galaxy for tool or workflow tests or runs, "
            "check that the history is in an 'ok' state before beginning tool "
            "or workflow execution."
        ),
    )


def required_tool_arg(allow_uris=False):
    """Decorate click method as requiring the path to a single tool."""
    arg_type_class = click.Path if not allow_uris else UriLike
    arg_type = arg_type_class(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
    if allow_uris:
        name = "uri"
        metavar = "TOOL_URI"
    else:
        name = "path"
        metavar = "TOOL_PATH"
    return click.argument(name, metavar=metavar, type=arg_type)


def paste_test_data_paths_option():
    return planemo_option(
        "--paste_test_data_paths/--no_paste_test_data_paths",
        is_flag=True,
        default=None,
        help=(
            "By default Planemo will use or not use Galaxy's path paste option to load "
            "test data into a history based on the engine type it is targeting. This can "
            "override the logic to explicitly enable or disable path pasting."
        ),
    )


def shed_install_option():
    return planemo_option(
        "--shed_install/--no_shed_install",
        is_flag=True,
        default=True,
        help=(
            "By default Planemo will attempt to install repositories needed for workflow "
            "testing. This may not be appropriate for production servers and so this can "
            "disabled by calling planemo with --no_shed_install."
        ),
    )


def install_tool_dependencies_option():
    return planemo_option(
        "--install_tool_dependencies/--no_install_tool_dependencies",
        is_flag=True,
        default=False,
        help=("Turn on installation of tool dependencies using classic toolshed packages."),
    )


def install_resolver_dependencies_option():
    return planemo_option(
        "--install_resolver_dependencies/--no_install_resolver_dependencies",
        is_flag=True,
        default=True,
        help=("Skip installing tool dependencies through resolver (e.g. conda)."),
    )


def install_repository_dependencies_option():
    return planemo_option(
        "--install_repository_dependencies/--no_install_repository_dependencies",
        is_flag=True,
        default=True,
        help=("Skip installing the repository dependencies."),
    )


def single_user_mode_option():
    return planemo_option(
        "galaxy_single_user",
        "--galaxy_single_user/--no_galaxy_single_user",
        is_flag=True,
        default=True,
        help=(
            "By default Planemo will configure Galaxy to run in single-user mode where there "
            "is just one user and this user is automatically logged it. Use --no_galaxy_single_user "
            "to prevent Galaxy from running this way."
        ),
    )


def required_workflow_arg():
    return click.argument(
        "workflow_identifier",
        metavar="WORKFLOW_PATH_OR_ID",
        type=str,
    )


def required_invocation_id_arg():
    return click.argument(
        "invocation_id",
        metavar="INVOCATION_ID",
        type=str,
    )


def invocation_export_format_arg():
    return click.argument(
        "export_format",
        default="rocrate.zip",
        metavar="model store format",
        type=str,
    )


def split_job_and_test():
    return click.option(
        "--split_test/--no_split_test", default=False, help="Write workflow job and test definitions to separate files."
    )


def from_invocation():
    return planemo_option(
        "--from_invocation/--from_uri",
        is_flag=True,
        default=False,
        help="Build a workflow test or job description from an invocation ID run on an external Galaxy."
        "A Galaxy URL and API key must also be specified. This allows test data to be downloaded"
        "and inputs and parameters defined automatically. Alternatively, the default is to build the"
        "descriptions from a provided workflow URI.",
    )


def required_job_arg():
    """Decorate click method as requiring the path to a single tool."""
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=False,
    )
    return click.argument("job_path", metavar="JOB_PATH", type=arg_type)


def required_runnable_arg():
    return click.argument(
        "runnable_identifier",
        metavar="RUNNABLE_PATH_OR_ID",
        type=str,
    )


def required_new_job_arg():
    arg_type = click.Path()
    return click.argument("new_job_path", metavar="NEW_JOB_PATH", type=arg_type)


def _optional_tools_default(ctx, param, value):
    if param.name in ["paths", "uris"] and len(value) == 0:
        return [os.path.abspath(os.getcwd())]
    else:
        return value


def optional_tools_or_packages_arg(multiple=False):
    """Decorate click method as optionally taking in the path to a tool
    or directory of tools or a Conda package. If no such argument is given
    the current working directory will be treated as a directory of tools.
    """
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="TARGET",
        nargs=nargs,
    )


class UriLike(click.Path):
    def convert(self, value, param, ctx):
        if "://" in value:
            return value
        else:
            return super().convert(value, param, ctx)


def optional_tools_arg(multiple=False, allow_uris=False, metavar="TOOL_PATH"):
    """Decorate click method as optionally taking in the path to a tool
    or directory of tools. If no such argument is given the current working
    directory will be treated as a directory of tools.
    """
    arg_type_class = click.Path if not allow_uris else UriLike
    arg_type = arg_type_class(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    )
    if allow_uris:
        name = "uris" if multiple else "uri"
    else:
        name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar=metavar,
        type=arg_type,
        nargs=nargs,
        callback=_optional_tools_default,
    )


class ProjectOrRepositry(click.Path):
    def __init__(self, **kwds):
        super().__init__(**kwds)

    def convert(self, value, param, ctx):
        if value and value.startswith("git:") or value.startswith("git+"):
            return value
        else:
            return super().convert(value, param, ctx)


def shed_project_arg(multiple=True):
    arg_type = ProjectOrRepositry(
        exists=True,
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
    )
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="PROJECT",
        type=arg_type,
        nargs=nargs,
        callback=_optional_tools_default,
    )


def recipe_arg(multiple=True):
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="RECIPE_DIR",
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=True,
            resolve_path=True,
        ),
        nargs=nargs,
        callback=_optional_tools_default,
    )


def optional_project_arg(exists=True, default="."):
    arg_type = click.Path(
        exists=exists,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    )
    return click.argument("path", metavar="PROJECT", default=default, type=arg_type)


def no_cleanup_option():
    return planemo_option("--no_cleanup", is_flag=True, help=("Do not cleanup temp files created for and by Galaxy."))


def docker_enable_option():
    return planemo_option("--docker/--no_docker", default=False, help=("Run Galaxy tools in Docker if enabled."))


def singularity_enable_option():
    return planemo_option(
        "--singularity/--no_singularity", default=False, help=("Run Galaxy tools in Singularity if enabled.")
    )


def docker_cmd_option():
    return planemo_option(
        "--docker_cmd",
        default=docker_util.DEFAULT_DOCKER_COMMAND,
        help="Command used to launch docker (defaults to docker).",
    )


def singularity_cmd_option():
    return planemo_option(
        "--singularity_cmd",
        default=singularity_util.DEFAULT_SINGULARITY_COMMAND,
        help="Command used to execute singularity (defaults to 'singularity').",
    )


def docker_sudo_option():
    return planemo_option("--docker_sudo/--no_docker_sudo", is_flag=True, help="Flag to use sudo when running docker.")


def singularity_sudo_option():
    return planemo_option(
        "--singularity_sudo/--no_singularity_sudo", is_flag=True, help="Flag to use sudo when running docker."
    )


def docker_sudo_cmd_option():
    return planemo_option(
        "--docker_sudo_cmd",
        help="sudo command to use when --docker_sudo is enabled " + "(defaults to sudo).",
        default=docker_util.DEFAULT_SUDO_COMMAND,
        use_global_config=True,
    )


def singularity_sudo_cmd_option():
    return planemo_option(
        "--singularity_sudo_cmd",
        help="sudo command to use when --singularity_sudo is enabled " + "(defaults to sudo).",
        default=singularity_util.DEFAULT_SUDO_COMMAND,
        use_global_config=True,
    )


def docker_host_option():
    return planemo_option(
        "--docker_host",
        help="Docker host to target when executing docker commands " + "(defaults to localhost).",
        use_global_config=True,
        default=docker_util.DEFAULT_HOST,
    )


def docker_run_extra_arguments_option():
    return planemo_option(
        "--docker_run_extra_arguments",
        help="Extra arguments to pass to docker run.",
        use_global_config=True,
        default="",
    )


def docker_config_options():
    return _compose(
        docker_cmd_option(),
        docker_sudo_option(),
        docker_host_option(),
        docker_sudo_cmd_option(),
        docker_run_extra_arguments_option(),
    )


def singularity_config_options():
    return _compose(
        singularity_cmd_option(),
        singularity_sudo_option(),
        singularity_sudo_cmd_option(),
    )


def galaxy_docker_options():
    return _compose(
        docker_enable_option(),
        docker_config_options(),
    )


def shed_owner_option():
    return planemo_option("--owner", help="Tool Shed repository owner (username).")


def shed_name_option():
    return planemo_option("--name", help="Tool Shed repository name (defaults to the inferred tool directory name).")


def validate_shed_target_callback(ctx, param, value):
    if value is None:
        ctx.fail("default_shed_target set to None, must specify a value for --shed_target to run this command.")
    return value


def shed_target_option():
    return planemo_option(
        "-t",
        "--shed_target",
        help="Tool Shed to target (this can be 'toolshed', 'testtoolshed', "
        "'local' (alias for http://localhost:9009/), an arbitrary url "
        "or mappings defined ~/.planemo.yml.",
        default=None,
        use_global_config=True,
        callback=validate_shed_target_callback,
    )


def shed_key_option():
    return planemo_option(
        "--shed_key",
        help=(
            "API key for Tool Shed access. An API key is required unless "
            "e-mail and password is specified. This key can be specified "
            "with either --shed_key or --shed_key_from_env."
        ),
    )


def shed_key_from_env_option():
    return planemo_option("--shed_key_from_env", help="Environment variable to read API key for Tool Shed access from.")


def shed_email_option():
    return planemo_option("--shed_email", help="E-mail for Tool Shed auth (required unless shed_key is specified).")


def shed_password_option():
    return planemo_option(
        "--shed_password", help="Password for Tool Shed auth (required unless shed_key is specified)."
    )


def shed_skip_upload():
    return planemo_option(
        "--skip_upload", is_flag=True, help=("Skip upload contents as part of operation, only update metadata.")
    )


def shed_skip_metadata():
    return planemo_option(
        "--skip_metadata",
        is_flag=True,
        help=("Skip metadata update as part of operation, only upload new contents."),
    )


def shed_message_option():
    return planemo_option("-m", "--message", help="Commit message for tool shed upload.")


def shed_force_create_option():
    return planemo_option(
        "--force_repository_creation",
        help=(
            "If a repository cannot be found for the specified user/repo "
            "name pair, then automatically create the repository in the "
            "toolshed."
        ),
        is_flag=True,
        default=False,
    )


def shed_check_diff_option():
    return planemo_option(
        "--check_diff",
        is_flag=True,
        help=(
            "Skip uploading if the shed_diff detects there would be no "
            "'difference' (only attributes populated by the shed would "
            "be updated.)"
        ),
    )


def shed_upload_options():
    return _compose(
        shed_message_option(),
        shed_force_create_option(),
        shed_check_diff_option(),
    )


def shed_realization_options():
    return _compose(
        shed_project_arg(multiple=True),
        recursive_shed_option(),
        shed_fail_fast_option(),
    )


def shed_repo_options():
    return _compose(
        shed_owner_option(),
        shed_name_option(),
    )


def shed_publish_options():
    """Common options for commands that require publishing to a
    a shed.
    """
    return _compose(
        shed_realization_options(),
        shed_repo_options(),
        shed_target_options(),
    )


def shed_read_options():
    """Common options that require read access to mapped repositories
    in a shed.
    """
    return _compose(
        shed_realization_options(),
        shed_repo_options(),
        shed_target_options(),
    )


def shed_target_options():
    """Common options for commands that require read-only
    interactions with a shed.
    """
    return _compose(
        shed_email_option(),
        shed_key_option(),
        shed_key_from_env_option(),
        shed_password_option(),
        shed_target_option(),
    )


def conda_target_options(include_local=True):
    return _compose(
        conda_prefix_option(),
        conda_exec_option(),
        conda_ensure_channels_option(),
        conda_use_local_option(),
    )


def github_namespace():
    return planemo_option(
        "--namespace",
        use_env_var=True,
        help=(
            "Organization or username under which to create or update workflow repository. "
            "Must be a valid github username or organization"
        ),
        default="iwc-workflows",
    )


def dry_run():
    return planemo_option(
        "--dry_run",
        help="Don't execute action, show preview of action.",
        is_flag=True,
        default=False,
    )


def github_branch():
    return planemo_option(
        "--github_branch",
        help="GitHub branch to use for the action. Default is main.",
        default="main",
    )


def galaxy_run_options():
    return _compose(
        galaxy_target_options(),
        galaxy_port_option(),
        galaxy_host_option(),
    )


def galaxy_config_options():
    return _compose(
        test_data_option(),
        tool_data_table_option(),
        dependency_resolvers_option(),
        brew_dependency_resolution(),
        shed_dependency_resolution(),
        no_dependency_resolution(),
        conda_target_options(),
        conda_dependency_resolution(),
        conda_auto_install_option(),
        conda_auto_init_option(),
        simultaneous_upload_option(),
        check_uploads_ok_option(),
        # Profile options...
        profile_option(),
        profile_database_options(),
        file_path_option(),
        database_connection_option(),
        postgres_database_storage_location_option(),
        shed_tools_conf_option(),
        shed_tools_directory_option(),
        single_user_mode_option(),
    )


def galaxy_target_options():
    return _compose(
        galaxy_root_option(),
        galaxy_python_version(),
        extra_tools_option(),
        install_galaxy_option(),
        galaxy_branch_option(),
        galaxy_source_option(),
        skip_venv_option(),
        no_cache_galaxy_option(),
        no_cleanup_option(),
        galaxy_email_option(),
        galaxy_docker_options(),
        mulled_containers_option(),
        galaxy_startup_timeout_option(),
        # Profile options...
        job_config_option(),
        tool_dependency_dir_option(),
        tool_data_path_option(),
    )


def pid_file_option():
    pid_file_type = click.Path(
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
    )
    return planemo_option(
        "--pid_file", type=pid_file_type, default=None, help="Location of pid file is executed with --daemon."
    )


def daemon_option():
    return planemo_option("--daemon", is_flag=True, help="Serve Galaxy process as a daemon.")


def profile_option(required=False):
    return planemo_option(
        "--profile",
        type=click.STRING,
        required=required,
        default=None,
        help=("Name of profile (created with the profile_create command) to use with this command."),
    )


def alias_option(required=False):
    return planemo_option("--alias", type=click.STRING, required=required, default=None, help=("Name of an alias."))


def galaxy_serve_options():
    return _compose(
        galaxy_run_options(),
        serve_engine_option(),
        non_strict_cwl_option(),
        docker_galaxy_image_option(),
        docker_extra_volume_option(),
        galaxy_config_options(),
        daemon_option(),
        pid_file_option(),
        ignore_dependency_problems_option(),
        install_prebuilt_client_option(),
        skip_client_build_option(),
        shed_install_option(),
    )


def training_topic_name_option():
    return planemo_option(
        "--topic_name",
        required=True,
        help="Name (directory name) of the topic to create or in which a tutorial should be created or updates",
    )


def training_topic_option():
    return _compose(
        training_topic_name_option(),
        planemo_option("--topic_title", default="Title of the topic", help="Title of the topic to create"),
        planemo_option("--topic_summary", default="Summary of the topic", help="Summary of the topic"),
        planemo_option(
            "--topic_target",
            type=click.Choice(["use", "admin-dev", "instructors"]),
            default="use",
            help="Target audience for the topic",
        ),
    )


def training_tutorial_name_option():
    return planemo_option("--tutorial_name", help="Name (directory name) of the tutorial to create or to modify")


def training_tutorial_name_req_option():
    return planemo_option("--tutorial_name", required=True, help="Name (directory name) of the tutorial to modify")


def training_zenodo_option():
    return planemo_option("--zenodo_link", help="Zenodo URL with the input data")


def training_tutorial_worflow_option():
    return _compose(
        planemo_option(
            "--workflow",
            type=click.Path(file_okay=True, resolve_path=True),
            help="Workflow of the tutorial (locally)",
            default=None,
        ),
        planemo_option("--galaxy_url", help="URL of a Galaxy instance with the workflow"),
        planemo_option("--galaxy_api_key", help="API key on the Galaxy instance with the workflow"),
        planemo_option("--workflow_id", help="ID of the workflow on the Galaxy instance"),
    )


def training_tutorial_option():
    return _compose(
        training_tutorial_name_option(),
        planemo_option("--tutorial_title", default="Title of the tutorial", help="Title of the tutorial"),
        planemo_option("--hands_on", is_flag=True, default=True, help="Add hands-on for the new tutorial"),
        planemo_option("--slides", is_flag=True, default=False, help="Add slides for the new tutorial"),
        training_tutorial_worflow_option(),
        training_zenodo_option(),
    )


def training_init_options():
    return _compose(training_topic_option(), training_tutorial_option())


def training_fill_data_library_options():
    return _compose(
        training_topic_name_option(),
        training_tutorial_name_req_option(),
        training_zenodo_option(),
    )


def training_generate_tuto_from_wf_options():
    return _compose(
        training_topic_name_option(), training_tutorial_name_req_option(), training_tutorial_worflow_option()
    )


def shed_fail_fast_option():
    return planemo_option(
        "--fail_fast",
        is_flag=True,
        default=False,
        help="If multiple repositories are specified and an error occurs "
        "stop immediately instead of processing remaining repositories.",
    )


def lint_biocontainers_option():
    return planemo_option(
        "biocontainer",
        "--biocontainer",
        "--biocontainers",
        is_flag=True,
        default=False,
        help="Check best practice BioContainer namespaces for a container definition applicable for this tool.",
    )


def report_level_option():
    return planemo_option(
        "--report_level",
        type=click.Choice(["all", "warn", "error"]),
        default="all",
    )


def report_xunit():
    return planemo_option(
        "--report_xunit",
        type=click.Path(file_okay=True, resolve_path=True),
        help="Output an XUnit report, useful for CI testing",
        default=None,
    )


def skip_options():
    return _compose(
        skip_option(),
        skip_file_option(),
    )


def skip_option():
    return planemo_option(
        "-s",
        "--skip",
        default=None,
        help=(
            "Comma-separated list of lint tests to skip (e.g. passing "
            "--skip 'citations,xml_order' would skip linting of citations "
            "and best-practice XML ordering."
        ),
    )


def skip_file_option():
    return planemo_option(
        "--skip_file",
        type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
        multiple=True,
        help=("File containing a list of lint tests to skip"),
    )


def fail_level_option():
    return planemo_option("--fail_level", type=click.Choice(["warn", "error"]), default="warn")


def recursive_shed_option():
    return recursive_option(
        "Recursively perform command for nested repository directories.",
    )


def recursive_option(help="Recursively perform command for subdirectories."):
    return planemo_option(
        "-r",
        "--recursive",
        is_flag=True,
        help=help,
    )


def merge_test_json():
    target_path = click.Path(
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
    )
    return click.argument(
        "input_paths",
        metavar="INPUT_PATHS",
        type=target_path,
        nargs=-1,
    )


def tool_test_json(var="path"):
    target_path = click.Path(
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
    )
    return click.argument(
        var,
        metavar="FILE_PATH",
        type=target_path,
        default="tool_test_output.json",
    )


def engine_options():
    return _compose(
        run_engine_option(),
        non_strict_cwl_option(),
        cwltool_no_container_option(),
        docker_galaxy_image_option(),
        docker_extra_volume_option(),
        ignore_dependency_problems_option(),
        shed_install_option(),
        install_tool_dependencies_option(),
        install_resolver_dependencies_option(),
        install_repository_dependencies_option(),
        galaxy_url_option(),
        galaxy_admin_key_option(),
        galaxy_user_key_option(),
        history_name(),
        history_id(),
        no_wait_option(),
    )


def test_report_options():
    return _compose(
        planemo_option(
            "--test_output",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            default="tool_test_output.html",
            help=("Output test report (HTML - for humans) defaults to tool_test_output.html."),
        ),
        planemo_option(
            "--test_output_text",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (Basic text - for display in CI)"),
            default=None,
        ),
        planemo_option(
            "--test_output_markdown",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (Markdown style - for humans & computers)"),
            default=None,
        ),
        planemo_option(
            "--test_output_markdown_minimal",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (Minimal markdown style - jost the table)"),
            default=None,
        ),
        planemo_option(
            "--test_output_xunit",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (xunit style - for CI systems"),
            default=None,
        ),
        planemo_option(
            "--test_output_junit",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (jUnit style - for CI systems"),
            default=None,
        ),
        planemo_option(
            "--test_output_allure",
            type=click.Path(file_okay=False, resolve_path=True),
            use_global_config=True,
            help=("Output test allure2 framework resutls"),
            default=None,
        ),
    )


def profile_name_argument():
    return click.argument(
        "profile_name",
        metavar="PROFILE_NAME",
        type=str,
    )


def database_identifier_argument():
    return click.argument(
        "identifier",
        metavar="IDENTIFIER",
        type=str,
    )


def postgres_option_callback(ctx, param, value):
    if value:
        ctx.fail("The `--postgres` option is deprecated, use `--database_type postgres` instead.")
    return value


def postgres_datatype_type_option():
    return click.option(
        "--postgres",
        is_flag=True,
        hidden=True,
        callback=postgres_option_callback,
    )


def postgres_database_storage_location_option():
    return planemo_option(
        "--postgres-storage-location",
        type=str,
        help="storage path for postgres database, used for local singularity postgres.",
        default=None,
        use_global_config=True,
    )


def no_wait_option():
    return planemo_option(
        "--no_wait",
        is_flag=True,
        default=False,
        prompt=False,
        help="After invoking a job or workflow, do not wait for completion.",
    )


def database_type_option():
    return planemo_option(
        "--database_type",
        default="auto",
        type=click.Choice(
            [
                "postgres",
                "postgres_docker",
                "postgres_singularity",
                "sqlite",
                "auto",
            ]
        ),
        use_global_config=True,
        help=(
            "Type of database to use for profile - "
            "'auto', 'sqlite', 'postgres', 'postgres_docker' , and postgres_singularity are available options. "
            "Use postgres to use an existing postgres server you user can "
            "access without a password via the psql command. Use postgres_docker "
            "to have Planemo manage a docker container running postgres. . Use "
            " postgres_singularity to have Planemo run postgres using singularity/apptainer. "
            "Data with postgres_docker is not yet persisted past when you restart "
            "the docker container launched by Planemo so be careful with this option."
        ),
    )


def database_source_options():
    """Database connection options for commands that utilize a database."""
    return _compose(
        planemo_option(
            "--postgres_psql_path",
            default="psql",
            use_global_config=True,
            help=("Name or or path to postgres client binary (psql)."),
        ),
        planemo_option(
            "--postgres_database_user",
            default="postgres",
            use_global_config=True,
            help=("Postgres username for managed development databases."),
        ),
        planemo_option(
            "--postgres_database_host",
            default=None,
            use_global_config=True,
            help=("Postgres host name for managed development databases."),
        ),
        planemo_option(
            "--postgres_database_port",
            default=None,
            use_global_config=True,
            help=("Postgres port for managed development databases."),
        ),
    )


def profile_database_options():
    return _compose(
        postgres_datatype_type_option(),
        database_type_option(),
        database_source_options(),
    )


def test_index_option():
    return planemo_option("--test_index", default=1, type=int, help="Select which test to check. Counting starts at 1")


def fail_fast_option():
    return planemo_option("--fail_fast", is_flag=True, help="Stop on first job failure.")


def test_output_options():
    return _compose(
        planemo_option(
            "--update_test_data",
            is_flag=True,
            help="Update test-data directory with job outputs (normally"
            " written to directory --job_output_files if specified.)",
        ),
        test_report_options(),
        planemo_option(
            "--test_output_json",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (planemo json) defaults to tool_test_output.json."),
            default="tool_test_output.json",
        ),
        planemo_option(
            "--job_output_files",
            type=click.Path(file_okay=False, resolve_path=True),
            help="Write job outputs to specified directory.",
            default=None,
        ),
        planemo_option(
            "--summary",
            type=click.Choice(["none", "minimal", "compact"]),
            default="minimal",
            help=(
                "Summary style printed to planemo's standard output (see "
                "output reports for more complete summary). Set to 'none' "
                "to disable completely."
            ),
        ),
        planemo_option(
            "--test_timeout",
            type=int,
            help="Maximum runtime of a single test in seconds.",
            default=DEFAULT_TOOL_TEST_WAIT,
        ),
    )


def test_options():
    return _compose(paste_test_data_paths_option(), test_output_options(), fail_fast_option())


def _compose(*functions):
    def compose2(f, g):
        return lambda x: f(g(x))

    return functools.reduce(compose2, functions)


def dependencies_script_options():
    return _compose(
        planemo_option(
            "--download_cache",
            type=click.Path(file_okay=False, resolve_path=True),
            use_global_config=True,
            help=("Directory to cache downloaded files, default is $DOWNLOAD_CACHE"),
            default=None,
        ),
    )


def filter_exclude_option():
    return planemo_option(
        "--exclude",
        type=click.Path(resolve_path=False),
        multiple=True,
        help="Paths to exclude.",
    )


def filter_exclude_from_option():
    return planemo_option(
        "--exclude_from",
        type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
        multiple=True,
        help="File of paths to exclude.",
    )


def filter_changed_in_commit_option():
    return planemo_option(
        "--changed_in_commit_range",
        help="Exclude paths unchanged in git commit range.",
    )


def ci_chunk_count_option():
    return planemo_option(
        "--chunk_count",
        type=int,
        help="Split output into chunks of this many item and print --chunk such group.",
        default=1,
    )


def ci_group_tools_option():
    return planemo_option("--group_tools", is_flag=True, help="Group tools of the same repository on a single line.")


def ci_chunk_option():
    return planemo_option(
        "--chunk",
        type=int,
        help=("When output is split into --chunk_count groups, output the group 0-indexedby this option."),
        default=0,
    )


def ci_output_option():
    return planemo_option(
        "--output",
        help="File to output to, or - for standard output.",
        default="-",
    )


def ci_find_options():
    return _compose(
        filter_exclude_option(),
        filter_exclude_from_option(),
        filter_changed_in_commit_option(),
        ci_chunk_count_option(),
        ci_chunk_option(),
        ci_output_option(),
    )


def workflow_output_artifact():
    return planemo_option(
        "-o",
        "--output",
        default=None,
        type=click.Path(
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    )


def tool_init_autopygen_option(prompt=False):
    return planemo_option(
        "--autopygen",
        type=click.STRING,
        prompt=prompt,
        help="Option for automatic generation of tool file,"
        " from python source code that uses argparse. "
        "Parameter is a path to source file containing definition of the parser",
    )


def tool_init_id_option(prompt=True):
    return planemo_option(
        "-i",
        "--id",
        type=click.STRING,
        prompt=prompt,
        help="Short identifier for new tool (no whitespace)",
    )


def tool_init_tool_option():
    return planemo_option(
        "-t",
        "--tool",
        default=None,
        type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True, resolve_path=True),
        help="Output path for new tool (default is <id>.xml)",
    )


def tool_init_name_option(prompt=True, help="Name for new tool (user facing)"):
    return planemo_option(
        "-n",
        "--name",
        type=click.STRING,
        prompt=prompt,
        help=help,
    )


def tool_init_version_option():
    return planemo_option(
        "--version",
        default="0.1.0",
        type=click.STRING,
        help="Tool XML version.",
    )


def tool_init_description_option():
    return planemo_option(
        "-d",
        "--description",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Short description for new tool (user facing)",
    )


def tool_init_command_option():
    return planemo_option(
        "-c",
        "--command",
        type=click.STRING,
        default=None,
        prompt=False,
        help=("Command potentially including cheetah variables ()(e.g. 'seqtk seq -A $input > $output')"),
    )


def tool_init_doi_option():
    return planemo_option(
        "--doi",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help=("Supply a DOI (http://www.doi.org/) easing citation of the tool for Galxy users (e.g. 10.1101/014043)."),
    )


def tool_init_test_case_option():
    return planemo_option(
        "--test_case",
        is_flag=True,
        default=None,
        prompt=False,
        help=("For use with --example_commmand, generate a tool test case from the supplied example."),
    )


def tool_init_macros_option():
    return planemo_option(
        "--macros",
        is_flag=True,
        default=None,
        prompt=False,
        help="Generate a macros.xml for reuse across many tools.",
    )


def tool_init_cite_url_option():
    return planemo_option(
        "--cite_url", type=click.STRING, default=None, multiple=True, prompt=False, help=("Supply a URL for citation.")
    )


def tool_init_input_option():
    return planemo_option(
        "--input",
        type=click.STRING,
        default=None,
        prompt=False,
        multiple=True,
        help="An input description (e.g. input.fasta)",
    )


def tool_init_output_option():
    return planemo_option(
        "--output",
        type=click.STRING,
        multiple=True,
        default=None,
        prompt=False,
        help=("An output location (e.g. output.bam), the Galaxy datatype is inferred from the extension."),
    )


def tool_init_help_text_option():
    return planemo_option(
        "--help_text",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Help text (reStructuredText)",
    )


def tool_init_help_from_command_option():
    return planemo_option(
        "--help_from_command",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Auto populate help from supplied command.",
    )


def tool_init_example_input_option():
    return planemo_option(
        "--example_input",
        type=click.STRING,
        default=None,
        prompt=False,
        multiple=True,
        help=("For use with --example_command, replace input file (e.g. 2.fastq with a data input parameter)."),
    )


def tool_init_example_output_option():
    return planemo_option(
        "--example_output",
        type=click.STRING,
        default=None,
        prompt=False,
        multiple=True,
        help=("For use with --example_command, replace input file (e.g. 2.fastq with a tool output)."),
    )


def tool_init_named_output_option():
    return planemo_option(
        "--named_output",
        type=click.STRING,
        multiple=True,
        default=None,
        prompt=False,
        help=(
            "Create a named output for use with command block for example "
            "specify --named_output=output1.bam and then use '-o $output1' "
            "in your command block."
        ),
    )


def tool_init_version_command_option():
    return planemo_option(
        "--version_command",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Command to print version (e.g. 'seqtk --version')",
    )


REQUIREMENT_HELP = "Add a tool requirement package (e.g. 'seqtk' or 'seqtk@1.68')."


def tool_init_requirement_option(help=REQUIREMENT_HELP):
    return planemo_option(
        "--requirement",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help=help,
    )


def tool_init_container_option():
    return planemo_option(
        "--container",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help="Add a Docker image identifier for this tool.",
    )


EXAMPLE_COMMAND_HELP = (
    "Example to command with paths to build Cheetah template from "
    "(e.g. 'seqtk seq -A 2.fastq > 2.fasta'). Option cannot be used "
    "with --command, should be used --example_input and "
    "--example_output."
)


def tool_init_example_command_option(help=EXAMPLE_COMMAND_HELP):
    return planemo_option(
        "--example_command",
        type=click.STRING,
        default=None,
        prompt=False,
        help=help,
    )


def mulled_conda_option():
    return planemo_option(
        "--mulled_conda_version",
        type=click.STRING,
        default=None,
        help=(
            "Install a specific version of Conda before running the command, by "
            "default the version that comes with the continuumio miniforge image "
            "will be used under Linux and under Mac OS X Conda will be upgraded to "
            "to work around a bug in 4.2."
        ),
    )


def mulled_namespace_option():
    return planemo_option(
        "--mulled_namespace",
        type=click.STRING,
        default="biocontainers",
        help=(
            "Build a mulled image with the specified namespace - defaults to "
            "biocontainers. Galaxy currently only recognizes images with the "
            "namespace biocontainers."
        ),
    )


def mulled_action_option():
    return planemo_option(
        "--mulled_command",
        type=click.STRING,
        default="build-and-test",
        help=("Mulled action to perform for targets - this defaults to 'build-and-test'."),
    )


def invocation_target_options():
    return _compose(
        required_invocation_id_arg(),
        galaxy_url_option(required=True),
        galaxy_user_key_option(required=True),
    )


def mulled_options():
    return _compose(
        mulled_conda_option(),
        mulled_namespace_option(),
        mulled_action_option(),
    )


def job_config_init_options():
    return _compose(
        docker_enable_option(),
        docker_config_options(),
        singularity_enable_option(),
        singularity_config_options(),
        extra_tools_option(),
        test_data_option(),
        tpv_option(),
        runner_target_option(),
        galaxy_version_option(),
    )
