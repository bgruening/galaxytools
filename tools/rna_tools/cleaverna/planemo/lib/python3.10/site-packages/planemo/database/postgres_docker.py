"""Module describes a :class:`DatabaseSource` for managed, dockerized postgres databases."""

import time

from galaxy.tool_util.deps import (
    docker_util,
    dockerfiles,
)
from galaxy.util import unicodify
from galaxy.util.commands import execute

from .interface import DatabaseSource
from .postgres import (
    _CommandBuilder,
    ExecutesPostgresSqlMixin,
)

DEFAULT_CONTAINER_NAME = "planemopostgres"
DEFAULT_POSTGRES_PASSWORD = "mysecretpassword"
DEFAULT_POSTGRES_PORT_EXPOSE = 15432


def docker_ps(args, **kwds):
    return docker_util.command_list("ps", args, **kwds)


def docker_exec(name, commands=[], **kwds):
    return docker_util.command_list("exec", [name] + commands, **kwds)


def is_running_container(name=DEFAULT_CONTAINER_NAME, **kwds):
    ps_command = docker_ps(["--format", "{{.Names}}"], **kwds)
    running_containers = unicodify(execute(ps_command))
    containers = running_containers.splitlines()
    return name in containers


def start_postgres_docker(
    name=DEFAULT_CONTAINER_NAME, password=DEFAULT_POSTGRES_PASSWORD, port=DEFAULT_POSTGRES_PORT_EXPOSE, **kwds
):
    run_command = docker_util.command_list(
        "run",
        ["-p", "%d:5432" % port, "--name", name, "-e", "POSTGRES_PASSWORD=%s" % password, "--rm", "-d", "postgres"],
        **kwds,
    )
    execute(run_command)


def stop_postgres_docker(name=DEFAULT_CONTAINER_NAME, **kwds):
    stop_command = docker_util.command_list("stop", [name], **kwds)
    execute(stop_command)


class DockerPostgresDatabaseSource(ExecutesPostgresSqlMixin, DatabaseSource):
    """Postgres database running inside a Docker container."""

    def __init__(self, **kwds):
        """Construct a postgres database source from planemo configuration."""
        self.psql_path = "psql"
        self.database_user = "postgres"
        self.database_password = DEFAULT_POSTGRES_PASSWORD
        self.database_host = "localhost"  # TODO: Make docker host
        self.database_port = DEFAULT_POSTGRES_PORT_EXPOSE
        self._kwds = kwds
        self._docker_host_kwds = dockerfiles.docker_host_args(**kwds)
        if not is_running_container(**self._docker_host_kwds):
            start_postgres_docker(**self._docker_host_kwds)
            # Hack to give docker a bit of time to boot up and allow psql to start.
            time.sleep(30)

    def sqlalchemy_url(self, identifier):
        """Return URL or form postgresql://username:password@localhost/mydatabase."""
        return "postgresql://%s:%s@%s:%d/%s" % (
            self.database_user,
            self.database_password,
            self.database_host,
            self.database_port,
            identifier,
        )

    def _psql_command_builder(self, *args):
        base_command = docker_exec(DEFAULT_CONTAINER_NAME, [self.psql_path], **self._docker_host_kwds)
        command_builder = _CommandBuilder(*base_command)
        # Print only tuples so output is easier to parse
        command_builder.append_command("--tuples-only")
        command_builder.append_command("--username", self.database_user)
        command_builder.append_command("-P", "pager=off")
        command_builder.extend_command(args)
        return command_builder


__all__ = ("DockerPostgresDatabaseSource",)
