"""Abstractions for setting up a Galaxy instance."""

import abc
import contextlib
import hashlib
import importlib.util
import json
import os
import random
import shlex
import shutil
import threading
import time
from string import Template
from tempfile import (
    mkdtemp,
    NamedTemporaryFile,
)
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Set,
    TYPE_CHECKING,
)

from galaxy.tool_util.deps import docker_util
from galaxy.util.commands import argv_to_str
from gxjobconfinit.generate import (
    build_job_config,
    ConfigArgs,
    create_docker_volumes,
    DevelopmentContext,
)
from packaging.version import parse as parse_version

from planemo import git
from planemo.config import OptionSource
from planemo.database import postgres_singularity
from planemo.deps import ensure_dependency_resolvers_conf_configured
from planemo.docker import docker_host_args
from planemo.galaxy.workflows import (
    get_toolshed_url_for_tool_id,
    remote_runnable_to_workflow_id,
)
from planemo.io import (
    communicate,
    kill_pid_file,
    shell,
    shell_join,
    stop_gravity,
    untar_to,
    wait_on,
    warn,
    write_file,
)
from planemo.mulled import build_involucro_context
from planemo.runnable import RunnableType
from planemo.shed import tool_shed_url
from .api import (
    DEFAULT_ADMIN_API_KEY,
    gi,
    user_api_key,
)
from .distro_tools import DISTRO_TOOLS_ID_TO_PATH
from .run import setup_venv
from .workflows import (
    find_tool_ids,
    import_workflow,
    install_shed_repos,
    MAIN_TOOLSHED_URL,
)

if TYPE_CHECKING:
    from planemo.runnable import Runnable

NO_TEST_DATA_MESSAGE = (
    "planemo couldn't find a target test-data directory, you should likely "
    "create a test-data directory or pass an explicit path using --test_data."
)

WEB_SERVER_CONFIG_TEMPLATE = """
[server:${server_name}]
use = egg:Paste#http
port = ${port}
host = ${host}
use_threadpool = True
threadpool_kill_thread_limit = 10800
[app:main]
paste.app_factory = galaxy.web.buildapp:app_factory
static_dir = static/
"""

TOOL_CONF_TEMPLATE = """<toolbox>
  <tool file="data_source/upload.xml" />
  ${tool_definition}
</toolbox>
"""

SHED_TOOL_CONF_TEMPLATE = """<?xml version="1.0"?>
<toolbox tool_path="${shed_tool_path}">
</toolbox>
"""

SHED_DATA_MANAGER_CONF_TEMPLATE = """<?xml version="1.0"?>
<data_managers>
</data_managers>
"""

EMPTY_JOB_METRICS_TEMPLATE = """<?xml version="1.0"?>
<job_metrics>
</job_metrics>
"""

FILE_SOURCES_TEMPLATE = """
- type: posix
  root: '${test_data_dir}'
  id: test_data_dir
  label: Test data directory
  doc: Test data directory for the runnables being tested
"""

TOOL_SHEDS_CONF = """<tool_sheds>
  <tool_shed name="Target Shed" url="${shed_target_url}" />
</tool_sheds>
"""

REFGENIE_CONFIG_TEMPLATE = """
config_version: %s
genome_folder: '%s'
genome_servers: ['http://refgenomes.databio.org']
genomes: null
"""

EMPTY_TOOL_CONF_TEMPLATE = """<toolbox></toolbox>"""
GX_TEST_TOOL_PATH = "$GALAXY_FUNCTIONAL_TEST_TOOLS"

DEFAULT_GALAXY_BRANCH = "master"
DEFAULT_GALAXY_SOURCE = "https://github.com/galaxyproject/galaxy"
CWL_GALAXY_SOURCE = "https://github.com/common-workflow-language/galaxy"

DATABASE_LOCATION_TEMPLATE = "sqlite:///%s?isolation_level=IMMEDIATE"

COMMAND_STARTUP_COMMAND = "./scripts/common_startup.sh ${COMMON_STARTUP_ARGS}"

CLEANUP_IGNORE_ERRORS = True
DEFAULT_GALAXY_BRAND = "Configured by Planemo"
DEFAULT_TOOL_INSTALL_TIMEOUT = 60 * 60 * 1
UNINITIALIZED = object()


@contextlib.contextmanager
def galaxy_config(ctx, runnables, **kwds):
    """Set up a ``GalaxyConfig`` in an auto-cleaned context."""
    c = local_galaxy_config
    if kwds.get("dockerize", False):
        c = docker_galaxy_config
    elif kwds.get("external", False):
        c = external_galaxy_config
    log_thread = None
    e = threading.Event()
    try:
        with c(ctx, runnables, **kwds) as config:
            if kwds.get("daemon"):
                log_thread = threading.Thread(target=read_log, args=(ctx, config.log_file, e))
                log_thread.daemon = True
                log_thread.start()
            yield config
    finally:
        if log_thread:
            e.set()
            log_thread.join(1)


def read_log(ctx, log_path, e: threading.Event):
    log_fh = None
    try:
        while not e.is_set():
            if os.path.exists(log_path):
                if not log_fh:
                    # Open in append so we start at the end of the log file
                    log_fh = open(log_path, "a+")
                log_lines = log_fh.read()
                if log_lines:
                    ctx.log(log_lines.rstrip())
            e.wait(1)
    finally:
        if log_fh:
            log_lines = log_fh.read()
            if log_lines:
                ctx.log(log_lines.rstrip())
            log_fh.close()


@contextlib.contextmanager
def docker_galaxy_config(ctx, runnables, for_tests=False, **kwds):
    """Set up a ``GalaxyConfig`` for Docker container."""
    test_data_dir = _find_test_data(runnables, **kwds)

    with _config_directory(ctx, **kwds) as config_directory:

        def config_join(*args):
            return os.path.join(config_directory, *args)

        ensure_dependency_resolvers_conf_configured(ctx, kwds, os.path.join(config_directory, "resolvers_conf.xml"))
        galaxy_root = kwds.get("galaxy_root")
        _handle_refgenie_config(config_directory, galaxy_root, kwds)

        shed_tool_conf = "config/shed_tool_conf.xml"
        all_tool_paths = _all_tool_paths(runnables, galaxy_root, kwds.get("extra_tools"))

        tool_directories = set()  # Things to mount...
        for tool_path in all_tool_paths:
            directory = os.path.dirname(os.path.normpath(tool_path))
            if os.path.exists(directory):
                tool_directories.add(directory)

        volumes = []
        for tool_directory in tool_directories:
            volumes.append(tool_directory)

        empty_tool_conf = config_join("empty_tool_conf.xml")

        tool_conf = config_join("tool_conf.xml")

        shed_tool_path = kwds.get("shed_tool_path") or config_join("shed_tools")
        _ensure_directory(shed_tool_path)

        sheds_config_path = _configure_sheds_config_file(ctx, config_directory, runnables, **kwds)
        port = _get_port(kwds)
        properties = _shared_galaxy_properties(config_directory, kwds, for_tests=for_tests)
        _handle_container_resolution(ctx, kwds, properties)
        master_api_key = _get_master_api_key(kwds)

        template_args = dict(
            shed_tool_path=shed_tool_path,
            tool_conf=tool_conf,
        )
        tool_config_file = f"{tool_conf},{shed_tool_conf}"

        _write_tool_conf(ctx, all_tool_paths, tool_conf)
        write_file(empty_tool_conf, EMPTY_TOOL_CONF_TEMPLATE)

        properties.update(
            dict(
                tool_config_file=tool_config_file,
                tool_sheds_config_file=sheds_config_path,
                migrated_tools_config=empty_tool_conf,
            )
        )

        server_name = f"planemo{random.randint(0, 100000)}"

        # Value substitutions in Galaxy properties - for consistency with
        # non-Dockerized version.
        template_args = dict()
        env = _build_env_for_galaxy(properties, template_args)
        env["NONUSE"] = "nodejs,proftp,reports"
        if ctx.verbose:
            env["GALAXY_LOGGING"] = "full"

        # TODO: setup FTP upload dir and disable FTP server in container.

        docker_target_kwds = docker_host_args(**kwds)
        volumes.append(config_directory)
        export_directory = kwds.get("export_directory", None)
        if export_directory is not None:
            volumes.append(f"{export_directory}:/export:rw")

        # TODO: Allow this to real Docker volumes and allow multiple.
        extra_volumes = kwds.get("docker_extra_volume") or []
        volumes.extend(extra_volumes)
        docker_volumes = create_docker_volumes(volumes)
        yield DockerGalaxyConfig(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            docker_target_kwds=docker_target_kwds,
            volumes=docker_volumes,
            export_directory=export_directory,
            kwds=kwds,
        )


@contextlib.contextmanager
def local_galaxy_config(ctx, runnables, for_tests=False, **kwds):
    """Set up a ``GalaxyConfig`` in an auto-cleaned context."""

    test_data_dir = _find_test_data(runnables, **kwds)
    tool_data_tables = _find_tool_data_table(runnables, test_data_dir=test_data_dir, **kwds)
    data_manager_config_paths = [r.data_manager_conf_path for r in runnables if r.data_manager_conf_path]
    galaxy_root = _find_galaxy_root(ctx, **kwds)
    install_galaxy = kwds.get("install_galaxy", False)
    if galaxy_root is not None:
        if os.path.isdir(galaxy_root) and not os.listdir(galaxy_root):
            os.rmdir(galaxy_root)
        if os.path.isdir(galaxy_root) and install_galaxy:
            raise Exception(f"{galaxy_root} is an existing non-empty directory, cannot install Galaxy again")

    # Duplicate block in docker variant above.
    if kwds.get("mulled_containers", False):
        if not kwds.get("docker", False):
            if ctx.get_option_source("docker") != OptionSource.cli:
                kwds["docker"] = True
            else:
                raise Exception("Specified no docker and mulled containers together.")
        conda_default_options = ("conda_auto_init", "conda_auto_install")
        use_conda_options = ("dependency_resolution", "conda_use_local", "conda_prefix", "conda_exec")
        if not any(kwds.get(_) for _ in use_conda_options) and all(
            ctx.get_option_source(_) == OptionSource.default for _ in conda_default_options
        ):
            # If using mulled_containers and default conda options disable conda resolution
            kwds["no_dependency_resolution"] = kwds["no_conda_auto_init"] = True

    with _config_directory(ctx, **kwds) as config_directory:

        def config_join(*args):
            return os.path.join(config_directory, *args)

        install_env = {}
        if kwds.get("galaxy_skip_client_build", True):
            install_env["GALAXY_SKIP_CLIENT_BUILD"] = "1"
        elif kwds.get("galaxy_install_prebuilt_client", True):
            install_env["GALAXY_INSTALL_PREBUILT_CLIENT"] = "1"
        if galaxy_root is None:
            galaxy_root = config_join("galaxy-dev")
        if not os.path.isdir(galaxy_root):
            _install_galaxy(ctx, galaxy_root, install_env, kwds)

        server_name = "main"
        # Once we don't have to support earlier than 18.01 - try putting these files
        # somewhere better than with Galaxy.
        log_file = f"{server_name}.log"
        pid_file = f"{server_name}.pid"
        ensure_dependency_resolvers_conf_configured(ctx, kwds, os.path.join(config_directory, "resolvers_conf.xml"))
        all_tool_paths = _all_tool_paths(runnables, galaxy_root=galaxy_root, extra_tools=kwds.get("extra_tools"))
        kwds["all_in_one_handling"] = True
        _handle_job_config_file(config_directory, server_name, test_data_dir, all_tool_paths, kwds)
        _handle_file_sources(config_directory, test_data_dir, kwds)
        _handle_refgenie_config(config_directory, galaxy_root, kwds)
        file_path = kwds.get("file_path") or config_join("files")
        _ensure_directory(file_path)

        tool_dependency_dir = kwds.get("tool_dependency_dir") or config_join("deps")
        _ensure_directory(tool_dependency_dir)

        shed_tool_conf = kwds.get("shed_tool_conf") or config_join("shed_tools_conf.xml")
        empty_tool_conf = config_join("empty_tool_conf.xml")

        tool_conf = config_join("tool_conf.xml")

        shed_data_manager_config_file = config_join("shed_data_manager_conf.xml")

        shed_tool_path = kwds.get("shed_tool_path") or config_join("shed_tools")
        _ensure_directory(shed_tool_path)

        sheds_config_path = _configure_sheds_config_file(ctx, config_directory, runnables, **kwds)

        database_location = config_join("galaxy.sqlite")
        master_api_key = _get_master_api_key(kwds)
        dependency_dir = os.path.join(config_directory, "deps")
        _ensure_directory(shed_tool_path)
        port = _get_port(kwds)
        template_args = dict(
            port=port,
            host=kwds.get("host", "127.0.0.1"),
            server_name=server_name,
            temp_directory=config_directory,
            shed_tool_path=shed_tool_path,
            database_location=database_location,
            tool_conf=tool_conf,
            debug=kwds.get("debug", "true"),
            id_secret=kwds.get("id_secret", hashlib.md5(str(time.time()).encode("utf-8")).hexdigest()),
            log_level="DEBUG" if ctx.verbose else "INFO",
        )
        tool_config_file = f"{tool_conf},{shed_tool_conf}"
        # Setup both galaxy_email and older test user test@bx.psu.edu
        # as admins for command_line, etc...
        properties = _shared_galaxy_properties(config_directory, kwds, for_tests=for_tests)
        properties.update(
            dict(
                server_name="main",
                enable_celery_tasks="true",
                ftp_upload_dir_template="${ftp_upload_dir}",
                ftp_upload_purge="false",
                ftp_upload_dir=test_data_dir or os.path.abspath("."),
                ftp_upload_site="Test Data",
                check_upload_content="false",
                tool_dependency_dir=dependency_dir,
                file_path=file_path,
                new_file_path="${temp_directory}/tmp",
                tool_config_file=tool_config_file,
                tool_sheds_config_file=sheds_config_path,
                manage_dependency_relationships="false",
                job_working_directory="${temp_directory}/job_working_directory",
                template_cache_path="${temp_directory}/compiled_templates",
                citation_cache_type="file",
                citation_cache_data_dir="${temp_directory}/citations/data",
                citation_cache_lock_dir="${temp_directory}/citations/lock",
                database_auto_migrate="true",
                enable_beta_tool_formats="true",
                id_secret="${id_secret}",
                log_level="${log_level}",
                debug="${debug}",
                watch_tools="auto",
                default_job_shell="/bin/bash",  # For conda dependency resolution
                tool_data_table_config_path=",".join(tool_data_tables) if tool_data_tables else None,
                data_manager_config_file=",".join(data_manager_config_paths)
                or None,  # without 'or None' may raise IOError in galaxy (see #946)
                integrated_tool_panel_config=("${temp_directory}/integrated_tool_panel_conf.xml"),
                migrated_tools_config=empty_tool_conf,
                test_data_dir=test_data_dir,  # TODO: make gx respect this
                shed_data_manager_config_file=shed_data_manager_config_file,
                outputs_to_working_directory="true",  # this makes Galaxy's files dir RO for dockerized testing
                object_store_store_by="uuid",
            )
        )
        _handle_container_resolution(ctx, kwds, properties)
        properties["database_connection"] = _database_connection(database_location, **kwds)
        if kwds.get("mulled_containers", False):
            properties["mulled_channels"] = kwds.get("conda_ensure_channels", "")

        _handle_kwd_overrides(properties, kwds)

        # TODO: consider following property
        # watch_tool = False
        # datatypes_config_file = config/datatypes_conf.xml
        # welcome_url = /static/welcome.html
        # logo_url = /
        # sanitize_all_html = True
        # serve_xss_vulnerable_mimetypes = False
        # track_jobs_in_database = None
        # retry_job_output_collection = 0

        env = _build_env_for_galaxy(properties, template_args)
        env.update(install_env)
        env["GALAXY_DEVELOPMENT_ENVIRONMENT"] = "1"
        # Following are needed in 18.01 to prevent Galaxy from changing log and pid.
        # https://github.com/galaxyproject/planemo/issues/788
        env["GALAXY_LOG"] = log_file
        env["GALAXY_PID"] = pid_file
        write_galaxy_config(
            galaxy_root=galaxy_root,
            properties=properties,
            env=env,
            kwds=kwds,
            template_args=template_args,
            config_join=config_join,
        )

        _write_tool_conf(ctx, all_tool_paths, tool_conf)
        write_file(empty_tool_conf, EMPTY_TOOL_CONF_TEMPLATE)

        shed_tool_conf_contents = _sub(SHED_TOOL_CONF_TEMPLATE, template_args)
        # Write a new shed_tool_conf.xml if needed.
        write_file(shed_tool_conf, shed_tool_conf_contents, force=False)

        write_file(shed_data_manager_config_file, SHED_DATA_MANAGER_CONF_TEMPLATE)

        yield LocalGalaxyConfig(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            galaxy_root,
            kwds,
        )


def write_galaxy_config(galaxy_root, properties, env, kwds, template_args, config_join):
    if get_galaxy_major_version(galaxy_root) < parse_version("22.01"):
        # Legacy .ini setup
        env["GALAXY_CONFIG_FILE"] = config_join("galaxy.ini")
        web_config = _sub(WEB_SERVER_CONFIG_TEMPLATE, template_args)
        write_file(env["GALAXY_CONFIG_FILE"], web_config)
    else:
        env["GALAXY_CONFIG_FILE"] = config_join("galaxy.yml")
        env["GRAVITY_STATE_DIR"] = config_join("gravity")
        with NamedTemporaryFile(suffix=".sock", delete=True) as nt:
            env["SUPERVISORD_SOCKET"] = nt.name
        write_file(
            env["GALAXY_CONFIG_FILE"],
            json.dumps(
                {
                    "galaxy": properties,
                    "gravity": {
                        "galaxy_root": galaxy_root,
                        "gunicorn": {
                            "bind": f"{kwds.get('host', 'localhost')}:{template_args['port']}",
                            "preload": False,
                        },
                        "gx_it_proxy": {
                            "enable": False,
                        },
                    },
                }
            ),
        )


def _expand_paths(galaxy_root: Optional[str], extra_tools: List[str]) -> List[str]:
    """Replace $GALAXY_FUNCTION_TEST_TOOLS with actual path."""
    if galaxy_root:
        extra_tools = [
            path if path != GX_TEST_TOOL_PATH else os.path.join(galaxy_root, "test/functional/tools")
            for path in extra_tools
        ]
    return extra_tools


def get_galaxy_major_version(galaxy_root):
    spec = importlib.util.spec_from_file_location(
        "__galaxy_version", os.path.join(galaxy_root, "lib", "galaxy", "version.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return parse_version(module.VERSION_MAJOR)


def get_refgenie_config(galaxy_root, refgenie_dir):
    config_version = 0.4
    if galaxy_root:
        version_major = get_galaxy_major_version(galaxy_root=galaxy_root)
        if version_major < parse_version("21.09"):
            config_version = 0.3
    return REFGENIE_CONFIG_TEMPLATE % (config_version, refgenie_dir)


def get_all_tool_path_from_kwds(runnables: List["Runnable"], **kwds) -> Set[str]:
    galaxy_root = kwds.get("galaxy_root")
    extra_tools = kwds.get("extra_tools")
    return _all_tool_paths(runnables, galaxy_root, extra_tools)


def _all_tool_paths(
    runnables: List["Runnable"], galaxy_root: Optional[str] = None, extra_tools: Optional[List[str]] = None
) -> Set[str]:
    extra_tools = extra_tools or []
    all_tool_paths = {r.path for r in runnables if r.has_tools and not r.data_manager_conf_path}
    extra_tools = _expand_paths(galaxy_root, extra_tools=extra_tools)
    all_tool_paths.update(extra_tools)
    for tool_id in get_tool_ids_for_runnables(runnables):
        tool_paths = DISTRO_TOOLS_ID_TO_PATH.get(tool_id)
        if tool_paths:
            if isinstance(tool_paths, str):
                tool_paths = [tool_paths]
            all_tool_paths.update(tool_paths)

    return all_tool_paths


def get_workflow_runnables(runnables: List["Runnable"]) -> List["Runnable"]:
    return [r for r in runnables if r.type == RunnableType.galaxy_workflow and r.has_path]


def get_tool_ids_for_runnables(runnables) -> List[str]:
    tool_ids = []
    for r in get_workflow_runnables(runnables):
        tool_ids.extend(find_tool_ids(r.path))
    return list(dict.fromkeys(tool_ids))


def _shared_galaxy_properties(config_directory, kwds, for_tests):
    """Setup properties useful for local and Docker Galaxy instances.

    Most things related to paths, etc... are very different between Galaxy
    modalities and many taken care of internally to the container in that mode.
    But this method sets up API stuff, tool, and job stuff that can be shared.
    """
    master_api_key = _get_master_api_key(kwds)
    user_email = _user_email(kwds)
    properties = {
        "master_api_key": master_api_key,
        "admin_users": f"{user_email},test@bx.psu.edu",
        "expose_dataset_path": "True",
        "collect_outputs_from": "job_working_directory",
        "allow_path_paste": "True",
        "check_migrate_tools": "False",
        "use_cached_dependency_manager": str(kwds.get("conda_auto_install", False)),
        "brand": kwds.get("galaxy_brand", DEFAULT_GALAXY_BRAND),
        "strict_cwl_validation": str(not kwds.get("non_strict_cwl", False)),
    }
    if kwds.get("no_cleanup", False):
        properties["cleanup_job"] = "never"
    else:
        properties["cleanup_job"] = "always"

    if kwds.get("galaxy_single_user", True):
        properties["single_user"] = user_email

    if for_tests:
        empty_dir = os.path.join(config_directory, "empty")
        _ensure_directory(empty_dir)
        properties["tour_config_dir"] = empty_dir
        properties["interactive_environment_plugins_directory"] = empty_dir
        properties["visualization_plugins_directory"] = empty_dir
        properties["refgenie_config_file"] = kwds.get("refgenie_config_file", "")
    return properties


@contextlib.contextmanager
def external_galaxy_config(ctx, runnables, for_tests=False, **kwds):
    yield BaseGalaxyConfig(
        ctx=ctx,
        galaxy_url=kwds.get("galaxy_url", None),
        master_api_key=_get_master_api_key(kwds),
        user_api_key=kwds.get("galaxy_user_key", None),
        runnables=runnables,
        kwds=kwds,
    )


def _get_master_api_key(kwds):
    master_api_key = kwds.get("galaxy_admin_key") or DEFAULT_ADMIN_API_KEY
    return master_api_key


def _get_port(kwds):
    port = int(kwds.get("port", 9090))
    return port


def _user_email(kwds):
    user_email = kwds.get("galaxy_email")
    return user_email


@contextlib.contextmanager
def _config_directory(ctx, **kwds):
    config_directory = kwds.get("config_directory", None)
    created_config_directory = False
    if not config_directory:
        created_config_directory = True
        config_directory = os.path.realpath(mkdtemp())
        ctx.vlog(f"Created directory for Galaxy configuration [{config_directory}]")
    try:
        yield config_directory
    finally:
        cleanup = not kwds.get("no_cleanup", False)
        if created_config_directory and cleanup:
            shutil.rmtree(config_directory, ignore_errors=True)


class GalaxyInterface(metaclass=abc.ABCMeta):
    """Abstraction around a Galaxy instance.

    Description of a Galaxy instance and how to interact with it - this could
    potentially be a remote, already running instance or an instance Planemo manages
    to execute some task(s).
    """

    @abc.abstractproperty
    def gi(self):
        """Return an admin bioblend Galaxy instance for API interactions."""

    @abc.abstractproperty
    def user_gi(self):
        """Return a user-backed bioblend Galaxy instance for API interactions."""

    @abc.abstractmethod
    def install_repo(self, *args, **kwds):
        """Install specified tool shed repository."""

    @abc.abstractproperty
    def tool_shed_client(self):
        """Return a admin bioblend tool shed client."""

    @abc.abstractmethod
    def wait_for_all_installed(self):
        """Wait for all queued up repositories installs to complete."""

    @abc.abstractmethod
    def install_workflows(self):
        """Install all workflows configured with these planemo arguments."""

    @abc.abstractmethod
    def workflow_id(self, path):
        """Get installed workflow API ID for input path."""

    @abc.abstractproperty
    def version_major(self):
        """Return target Galaxy version."""

    @abc.abstractproperty
    def user_api_config(self):
        """Return the API indicated configuration for user session.

        Calling .config.get_config() with admin GI session would yield
        a different object (admins have different view of Galaxy's
        configuration).
        """

    @property
    def user_is_admin(self):
        return self.user_api_config["is_admin_user"]


class GalaxyConfig(GalaxyInterface, metaclass=abc.ABCMeta):
    """Specialization of GalaxyInterface for Galaxy instances Planemo manages itself.

    This assumes more than an API connection is available - Planemo needs to be able to
    start and stop the Galaxy instance, recover logs, etc... There are currently two
    implementations - a locally executed Galaxy and one running inside a Docker containe
    """

    @abc.abstractproperty
    def kill(self):
        """Stop the running instance."""

    @abc.abstractmethod
    def startup_command(self, ctx, **kwds):
        """Return a shell command used to startup this instance.

        Among other common planmo kwds, this should respect the
        ``daemon`` keyword.
        """

    @abc.abstractproperty
    def log_contents(self):
        """Retrieve text of log for running Galaxy instance."""

    @abc.abstractmethod
    def cleanup(self):
        """Cleanup allocated resources to run this instance."""

    @abc.abstractproperty
    def use_path_paste(self):
        """Use path paste to upload data.

        This will only be an option if the target user key is an
        admin user key.
        """


class BaseGalaxyConfig(GalaxyInterface):
    def __init__(
        self,
        ctx,
        galaxy_url,
        master_api_key,
        user_api_key,
        runnables,
        kwds,
    ):
        self._ctx = ctx
        self.galaxy_url = galaxy_url
        self.master_api_key = master_api_key
        self._user_api_key = user_api_key
        self.runnables = runnables
        self._kwds = kwds
        self._workflow_ids = {}
        self.installed_repos = {}
        self.updated_repos = {}

        self._target_version = UNINITIALIZED
        self._target_user_config = UNINITIALIZED

    @property
    def gi(self):
        assert self.galaxy_url
        return gi(url=self.galaxy_url, key=self.master_api_key)

    @property
    def user_gi(self):
        user_api_key = self.user_api_key
        assert user_api_key
        return self._gi_for_key(user_api_key)

    @property
    def user_api_key(self):
        # TODO: thread-safe
        if self._user_api_key is None:
            # TODO: respect --galaxy_email - seems like a real bug
            self._user_api_key = user_api_key(self.gi)

        return self._user_api_key

    def _gi_for_key(self, key):
        assert self.galaxy_url
        return gi(url=self.galaxy_url, key=key)

    def install_repo(self, *args, **kwds):
        self.tool_shed_client.install_repository_revision(*args, **kwds)

    @property
    def tool_shed_client(self):
        return self.gi.toolShed

    def wait_for_all_installed(self):
        def status_ready(repo):
            status = repo["status"]
            if status in ["Installing", "New"]:
                return None
            if status == "Installed":
                return True
            raise Exception(f"Error installing repo status is {status}")

        def ready():
            repos = self.tool_shed_client.get_repositories()
            ready = all(map(status_ready, repos))
            return ready or None

        wait_on(ready, "galaxy tool installation", timeout=DEFAULT_TOOL_INSTALL_TIMEOUT)

    def install_workflows(self):
        for runnable in self.runnables:
            if runnable.type.name in ["galaxy_workflow", "cwl_workflow"] and not runnable.is_remote_workflow_uri:
                self._install_workflow(runnable)

    def _install_workflow(self, runnable):
        if self._kwds.get("shed_install") and (
            self._kwds.get("engine") != "external_galaxy" or self._kwds.get("galaxy_admin_key")
        ):
            workflow_repos = install_shed_repos(
                runnable,
                self.gi,
                self._kwds.get("ignore_dependency_problems", False),
                self._kwds.get("install_tool_dependencies", False),
                self._kwds.get("install_resolver_dependencies", True),
                self._kwds.get("install_repository_dependencies", True),
                self._kwds.get("install_most_recent_revision", False),
            )
            self.installed_repos[runnable.path], self.updated_repos[runnable.path] = workflow_repos
        default_from_path = self._kwds.get("workflows_from_path", False)
        # TODO: Allow serialization so this doesn't need to assume a
        # shared filesystem with Galaxy server.
        from_path = default_from_path or (runnable.type.name == "cwl_workflow")
        workflow = import_workflow(runnable.path, admin_gi=self.gi, user_gi=self.user_gi, from_path=from_path)
        self._workflow_ids[runnable.path] = workflow["id"]

    def workflow_id_for_runnable(self, runnable):
        if runnable.is_remote_workflow_uri:
            workflow_id = remote_runnable_to_workflow_id(runnable)
        else:
            workflow_id = self.workflow_id(runnable.path)
        return workflow_id

    def workflow_id(self, path):
        return self._workflow_ids[path]

    @property
    def use_path_paste(self):
        option = self._kwds.get("paste_test_data_paths")
        if option is None:
            return self.default_use_path_paste
        else:
            return option

    @property
    def default_use_path_paste(self):
        return False

    @property
    def version_major(self) -> str:
        """Return target Galaxy version."""
        if self._target_version is UNINITIALIZED:
            self._target_version = self.user_gi.config.get_version()["version_major"]
        return self._target_version

    @property
    def user_api_config(self):
        """Return the API indicated configuration for user session."""
        if self._target_user_config is UNINITIALIZED:
            self._target_user_config = self.user_gi.config.get_config()
        return self._target_user_config


class BaseManagedGalaxyConfig(BaseGalaxyConfig):
    def __init__(
        self,
        ctx,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
        runnables,
        kwds,
    ):
        galaxy_url = f"http://localhost:{port}"
        super().__init__(
            ctx=ctx,
            galaxy_url=galaxy_url,
            master_api_key=master_api_key,
            user_api_key=None,
            runnables=runnables,
            kwds=kwds,
        )
        self.config_directory = config_directory
        self.env = env
        self.test_data_dir = test_data_dir
        self.port = port
        self.server_name = server_name

    @property
    def log_file(self):
        """Log file used when planemo serves this Galaxy instance."""
        file_name = f"{self.server_name}.log"
        return file_name


class DockerGalaxyConfig(BaseManagedGalaxyConfig):
    """A :class:`GalaxyConfig` description of a Dockerized Galaxy instance."""

    def __init__(
        self,
        ctx,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
        runnables,
        docker_target_kwds,
        volumes,
        export_directory,
        kwds,
    ):
        super().__init__(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            kwds,
        )
        self.docker_target_kwds = docker_target_kwds
        self.volumes = volumes
        self.export_directory = export_directory

    def kill(self):
        """Kill planemo container..."""
        kill_command = docker_util.kill_command(self.server_name, **self.docker_target_kwds)
        return shell(kill_command)

    def startup_command(self, ctx, **kwds):
        """Return a shell command used to startup this instance.

        Among other common planmo kwds, this should respect the
        ``daemon`` keyword.
        """
        daemon = kwds.get("daemon", False)
        daemon_str = "" if not daemon else " -d"
        docker_run_extras = f"-p {self.port}:80{daemon_str}"
        env_directives = [f"{k}='{v}'" for k, v in self.env.items()]
        image = kwds.get("docker_galaxy_image", "bgruening/galaxy-stable")
        run_command = docker_util.build_docker_run_command(
            "",
            image,
            interactive=False,
            env_directives=env_directives,
            working_directory=None,
            name=self.server_name,
            run_extra_arguments=docker_run_extras,
            set_user=False,
            volumes=self.volumes,
            **self.docker_target_kwds,
        )
        chmod_command = [
            "chmod",
            "-R",
            "o+rwx",
            self.config_directory,
        ]
        if self.export_directory:
            chmod_command.append(self.export_directory)

        return shell_join(
            argv_to_str(chmod_command),
            run_command,
        )

    @property
    def log_contents(self):
        logs_command = docker_util.logs_command(self.server_name, **self.docker_target_kwds)
        output, _ = communicate(logs_command)
        return output

    def cleanup(self):
        shutil.rmtree(self.config_directory, CLEANUP_IGNORE_ERRORS)


class LocalGalaxyConfig(BaseManagedGalaxyConfig):
    """A local, non-containerized implementation of :class:`GalaxyConfig`."""

    def __init__(
        self,
        ctx,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
        runnables,
        galaxy_root,
        kwds,
    ):
        super().__init__(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            kwds,
        )
        self.galaxy_root = galaxy_root
        self._virtual_env_locs = []

    @property
    def virtual_env_dir(self):
        loc = None
        for loc in self._virtual_env_locs:
            if not os.path.isabs(loc):
                loc = os.path.join(self.galaxy_root, loc)
            if os.path.isdir(loc):
                break
        return loc

    @property
    def gravity_state_dir(self):
        return self.env["GRAVITY_STATE_DIR"]

    def kill(self):
        if self._ctx.verbose:
            shell(["ps", "ax"])
            exists = os.path.exists(self.pid_file)
            print(f"Killing pid file [{self.pid_file}]")
            print(f"pid_file exists? [{exists}]")
            if exists:
                with open(self.pid_file) as f:
                    print(f"pid_file contents are [{f.read()}]")
        if self.env.get("GRAVITY_STATE_DIR"):
            stop_gravity(
                virtual_env=self.virtual_env_dir or os.path.join(self.galaxy_root, ".venv"),
                gravity_state_dir=self.gravity_state_dir,
                env=self.env,
            )
        kill_pid_file(self.pid_file)

    def startup_command(self, ctx, **kwds):
        """Return a shell command used to startup this instance.

        Among other common planemo kwds, this should respect the
        ``daemon`` keyword.
        """
        daemon = kwds.get("daemon", False)
        # TODO: Allow running dockerized Galaxy here instead.
        setup_venv_command = setup_venv(ctx, kwds, self)
        run_script = f"{shlex.quote(os.path.join(self.galaxy_root, 'run.sh'))} $COMMON_STARTUP_ARGS"
        if daemon:
            run_script += " --daemon"
        else:
            run_script += f" --server-name {shlex.quote(self.server_name)}"
        cd_to_galaxy_command = ["cd", self.galaxy_root]
        return shell_join(
            cd_to_galaxy_command,
            setup_venv_command,
            run_script,
        )

    @property
    def log_file(self):
        """Log file used when planemo serves this Galaxy instance."""
        file_name = f"{self.server_name}.log"
        return os.path.join(self.galaxy_root, file_name)

    @property
    def pid_file(self):
        pid_file_name = f"{self.server_name}.pid"
        return os.path.join(self.galaxy_root, pid_file_name)

    @property
    def log_contents(self):
        if not os.path.exists(self.log_file):
            return ""
        with open(self.log_file) as f:
            return f.read()

    def cleanup(self):
        shutil.rmtree(self.config_directory, CLEANUP_IGNORE_ERRORS)

    @property
    def default_use_path_paste(self):
        # If Planemo started a local, native Galaxy instance assume files URLs can be
        # pasted.
        return self.user_is_admin


def _database_connection(database_location, **kwds):
    if "database_type" in kwds and kwds["database_type"] == "postgres_singularity":
        default_connection = postgres_singularity.DEFAULT_CONNECTION_STRING
    else:
        default_connection = DATABASE_LOCATION_TEMPLATE % database_location
    database_connection = kwds.get("database_connection") or default_connection
    return database_connection


def _find_galaxy_root(ctx, **kwds):
    root_prop = "galaxy_root"
    cwl = kwds.get("cwl", False)
    if cwl:
        root_prop = "cwl_galaxy_root"
    galaxy_root = kwds.get(root_prop, None)
    if galaxy_root:
        return galaxy_root
    else:
        par_dir = os.getcwd()
        while True:
            run = os.path.join(par_dir, "run.sh")
            config = os.path.join(par_dir, "config")
            if os.path.isfile(run) and os.path.isdir(config):
                return par_dir
            new_par_dir = os.path.dirname(par_dir)
            if new_par_dir == par_dir:
                break
            par_dir = new_par_dir
    return None


def _find_test_data(runnables, **kwds):
    test_data_search_path = "."
    runnables = [r for r in runnables if r.has_tools]
    if len(runnables) > 0:
        test_data_search_path = runnables[0].test_data_search_path

    # Find test data directory associated with path.
    test_data = kwds.get("test_data", None)
    if test_data:
        return os.path.abspath(test_data)
    else:
        test_data = _search_tool_path_for(test_data_search_path, "test-data")
        if test_data:
            return test_data
    warn(NO_TEST_DATA_MESSAGE)
    return None


def _find_tool_data_table(runnables, test_data_dir, **kwds) -> Optional[List[str]]:
    tool_data_search_path = "."
    runnables = [r for r in runnables if r.has_tools]
    if len(runnables) > 0:
        tool_data_search_path = runnables[0].tool_data_search_path

    tool_data_table = kwds.get("tool_data_table", None)
    if tool_data_table:
        return [os.path.abspath(table_path) for table_path in tool_data_table]
    else:
        extra_paths = [test_data_dir] if test_data_dir else []
        tool_data_table = _search_tool_path_for(
            tool_data_search_path,
            "tool_data_table_conf.xml.test",
            extra_paths,
        ) or _search_tool_path_for(  # if all else fails just use sample
            tool_data_search_path, "tool_data_table_conf.xml.sample"
        )
        if tool_data_table:
            return [tool_data_table]
    return None


def _search_tool_path_for(path, target, extra_paths=None):
    """Check for presence of a target in different artifact directories."""
    if extra_paths is None:
        extra_paths = []
    if not os.path.isdir(path):
        tool_dir = os.path.dirname(path)
    else:
        tool_dir = path
    possible_dirs = [tool_dir, "."] + extra_paths
    for possible_dir in possible_dirs:
        possible_path = os.path.join(possible_dir, target)
        if os.path.exists(possible_path):
            return os.path.abspath(possible_path)
    return None


def get_tool_sheds_conf_for_runnables(runnables: Optional[List["Runnable"]]) -> Optional[str]:
    if runnables:
        tool_ids = get_tool_ids_for_runnables(runnables)
        return get_shed_tools_conf_string_for_tool_ids(tool_ids)
    return None


def get_shed_tools_conf_string_for_tool_ids(tool_ids: List[str]) -> str:
    tool_shed_urls = set(get_toolshed_url_for_tool_id(tool_id) for tool_id in tool_ids if tool_id)
    # always add main toolshed
    tool_shed_urls.add(MAIN_TOOLSHED_URL)
    cleaned_tool_shed_urls = set(_ for _ in tool_shed_urls if _ is not None)
    TOOL_SHEDS_CONF_TEMPLATE = Template("""<tool_sheds>${tool_shed_lines}</tool_sheds>""")
    tool_sheds: List[str] = []
    # sort tool_shed_urls from shortest to longest, as https://github.com/galaxyproject/galaxy/blob/c7cb47a1b18ccd5b39075a705bbd2f34572755fe/lib/galaxy/util/tool_shed/tool_shed_registry.py#L106-L118
    # has a bug where a toolshed that is an exact substring of another registered toolshed would wrongly be selected.
    for shed_url in sorted(cleaned_tool_shed_urls, key=lambda url: len(url)):
        tool_sheds.append(f'<tool_shed name="{shed_url.split("://")[-1]}" url="{shed_url}" />')
    return TOOL_SHEDS_CONF_TEMPLATE.substitute(tool_shed_lines="".join(tool_sheds))


def _configure_sheds_config_file(ctx, config_directory, runnables, **kwds):
    # Find tool sheds to add to config
    contents = get_tool_sheds_conf_for_runnables(runnables)
    if not contents:
        if "shed_target" not in kwds:
            kwds = kwds.copy()
            kwds["shed_target"] = "toolshed"
        shed_target_url = tool_shed_url(ctx, **kwds)
        contents = _sub(TOOL_SHEDS_CONF, {"shed_target_url": shed_target_url})
    tool_sheds_conf = os.path.join(config_directory, "tool_sheds_conf.xml")
    write_file(tool_sheds_conf, contents)
    return tool_sheds_conf


def _tool_conf_entry_for(tool_paths):
    tool_definitions = ""
    for tool_path in tool_paths:
        if os.path.isdir(tool_path):
            tool_definitions += f"""<tool_dir dir="{tool_path}" />"""
        else:
            tool_definitions += f"""<tool file="{tool_path}" />"""
    return tool_definitions


def _install_galaxy(ctx, galaxy_root, env, kwds):
    if not kwds.get("no_cache_galaxy", False):
        _install_galaxy_via_git(ctx, galaxy_root, env, kwds)
    else:
        _install_galaxy_via_download(ctx, galaxy_root, env, kwds)


def _install_galaxy_via_download(ctx, galaxy_root, env, kwds):
    branch = _galaxy_branch(kwds)
    source = _galaxy_source(kwds)
    if source.startswith("https://github.com/"):
        source = source[len("https://github.com/") :]
    untar_to(
        f"https://codeload.github.com/{source}/tar.gz/{branch}",
        tar_args=["--strip-components", "1", "-xvzf", "-", "galaxy-" + branch.replace("/", "-")],
        dest_dir=galaxy_root,
    )
    _install_with_command(ctx, galaxy_root, env, kwds)


def _install_galaxy_via_git(ctx, galaxy_root, env, kwds):
    gx_repo = _ensure_galaxy_repository_available(ctx, kwds)
    branch = _galaxy_branch(kwds)
    command = git.command_clone(ctx, gx_repo, galaxy_root, branch=branch, depth=1)
    exit_code = shell(command, env=env)
    if exit_code != 0:
        raise Exception("Failed to clone Galaxy via git")
    _install_with_command(ctx, galaxy_root, env, kwds)


def _galaxy_branch(kwds):
    branch = kwds.get("galaxy_branch", None)
    if branch is None:
        cwl = kwds.get("cwl", False)
        branch = "cwl-1.0" if cwl else None
    if branch is None:
        branch = DEFAULT_GALAXY_BRANCH

    return branch


def _galaxy_source(kwds):
    source = kwds.get("galaxy_source", None)
    if source is None:
        cwl = kwds.get("cwl", False)
        source = CWL_GALAXY_SOURCE if cwl else None
    if source is None:
        source = DEFAULT_GALAXY_SOURCE

    return source


def _install_with_command(ctx, galaxy_root, env, kwds):
    setup_venv_command = setup_venv(ctx, kwds)
    install_cmd = shell_join(
        setup_venv_command,
        COMMAND_STARTUP_COMMAND,
    )
    communicate(install_cmd, default_err_msg="Failed to install Galaxy via command", cwd=galaxy_root, env=env)
    if not os.path.exists(galaxy_root):
        raise Exception(f"Failed to create Galaxy directory [{galaxy_root}]")
    if not os.path.exists(os.path.join(galaxy_root, "lib")):
        raise Exception(f"Failed to create Galaxy directory [{galaxy_root}], lib missing")


def _ensure_galaxy_repository_available(ctx, kwds):
    workspace = ctx.workspace
    cwl = kwds.get("cwl", False)
    galaxy_source = kwds.get("galaxy_source")
    if galaxy_source and galaxy_source != DEFAULT_GALAXY_SOURCE:
        sanitized_repo_name = "".join(c if c.isalnum() else "_" for c in kwds["galaxy_source"]).rstrip()[:255]
        gx_repo = os.path.join(workspace, f"gx_repo_{sanitized_repo_name}")
    else:
        gx_repo = os.path.join(workspace, "gx_repo")
    if cwl:
        gx_repo += "_cwl"
    if os.path.exists(gx_repo):
        # Convert the git repository from bare to mirror, if needed
        shell(["git", "--git-dir", gx_repo, "config", "remote.origin.fetch", "+refs/*:refs/*"])
        shell(["git", "--git-dir", gx_repo, "config", "remote.origin.mirror", "true"])
        # Attempt remote update - but don't fail if not interweb, etc...
        shell(f"git --git-dir {gx_repo} remote update >/dev/null 2>&1")
    else:
        remote_repo = _galaxy_source(kwds)
        command = git.command_clone(ctx, remote_repo, gx_repo, mirror=True)
        shell(command)
    return gx_repo


def _build_env_for_galaxy(properties, template_args):
    env = {}
    for key, value in properties.items():
        if value is not None:  # Do not override None with empty string
            var = f"GALAXY_CONFIG_OVERRIDE_{key.upper()}"
            value = _sub(value, template_args)
            env[var] = value
    return env


def _handle_job_config_file(
    config_directory: str,
    server_name: str,
    test_data_dir: Optional[str],
    all_tool_paths: Set[str],
    kwds: Dict[str, Any],
):
    job_config_file = kwds.get("job_config_file", None)
    if not job_config_file:
        dev_context = DevelopmentContext(
            test_data_dir,
            list(all_tool_paths),
        )
        init_config = ConfigArgs.from_dict(**kwds)
        conf_contents = build_job_config(init_config, dev_context)
        job_config_file = os.path.join(
            config_directory,
            "job_conf.yml",
        )
        write_file(job_config_file, conf_contents)
    kwds["job_config_file"] = job_config_file


def _write_tool_conf(ctx, tool_paths, tool_conf_path):
    tool_definition = _tool_conf_entry_for(tool_paths)
    tool_conf_template_kwds = dict(tool_definition=tool_definition)
    tool_conf_contents = _sub(TOOL_CONF_TEMPLATE, tool_conf_template_kwds)
    write_file(tool_conf_path, tool_conf_contents)
    ctx.vlog(
        "Writing tool_conf to path %s with contents [%s]",
        tool_conf_path,
        tool_conf_contents,
    )


def _handle_container_resolution(ctx, kwds, galaxy_properties):
    if kwds.get("mulled_containers", False):
        galaxy_properties["enable_beta_mulled_containers"] = "True"
        involucro_context = build_involucro_context(ctx, **kwds)
        galaxy_properties["involucro_auto_init"] = "False"  # Use planemo's
        galaxy_properties["involucro_path"] = involucro_context.involucro_bin


def _handle_file_sources(config_directory, test_data_dir, kwds):
    file_sources_conf = os.path.join(config_directory, "file_sources_conf.yml")
    file_sources_conf_contents = _sub(FILE_SOURCES_TEMPLATE, {"test_data_dir": test_data_dir})
    write_file(file_sources_conf, file_sources_conf_contents)
    kwds["file_sources_config_file"] = file_sources_conf


def _handle_refgenie_config(config_directory, galaxy_root, kwds):
    refgenie_dir = os.path.join(config_directory, "refgenie")
    _ensure_directory(refgenie_dir)
    refgenie_config_file = os.path.join(refgenie_dir, "genome_config.yaml")
    refgenie_config = get_refgenie_config(galaxy_root=galaxy_root, refgenie_dir=refgenie_dir)
    with open(refgenie_config_file, "w") as fh:
        fh.write(refgenie_config)
    kwds["refgenie_config_file"] = refgenie_config_file


def _handle_kwd_overrides(properties, kwds):
    kwds_gx_properties = [
        "tool_data_path",
        "job_config_file",
        "job_metrics_config_file",
        "dependency_resolvers_config_file",
    ]
    for prop in kwds_gx_properties:
        val = kwds.get(prop, None)
        if val:
            properties[prop] = val


def _sub(template, args):
    if template is None:
        return ""
    return Template(template).safe_substitute(args)


def _ensure_directory(path):
    if path is not None and not os.path.exists(path):
        os.makedirs(path)


__all__ = (
    "DATABASE_LOCATION_TEMPLATE",
    "galaxy_config",
)
