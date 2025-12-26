"""Module describes a :class:`DatabaseSource` for local postgres databases."""

import subprocess

from galaxy.util import unicodify

from planemo.io import communicate
from .interface import DatabaseSource


class ExecutesPostgresSqlMixin:
    def list_databases(self):
        """Use `psql --list` to generate a list of identifiers."""
        command_builder = self._psql_command_builder("--list")
        stdout = unicodify(self._communicate(command_builder))
        output_lines = stdout.splitlines()
        identifiers = []
        for line in output_lines:
            identifiers.append(line.split("|")[0].strip())
        return [i for i in identifiers if i]

    def create_database(self, identifier):
        """Use `psql -c "create database"` to create a database."""
        sql = "create database %s;" % identifier
        self._run_sql_command(sql)

    def delete_database(self, identifier):
        """Use `psql -c "drop database"` to delete a database."""
        sql = "drop database %s;" % identifier
        self._run_sql_command(sql)

    def _run_sql_command(self, sql):
        # communicate is just joining commands so we need to modify the
        # sql as an argument - it shouldn't do this.
        sql_arg = "%s" % sql
        command_builder = self._psql_command_builder("--command", sql_arg)
        self._communicate(command_builder)

    def _communicate(self, command_builder):
        stdout, _ = communicate(
            command_builder.command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return stdout


class LocalPostgresDatabaseSource(ExecutesPostgresSqlMixin, DatabaseSource):
    """Local postgres database source managed through psql application."""

    def __init__(self, **kwds):
        """Construct a postgres database source from planemo configuration."""
        self.psql_path = kwds.get("postgres_psql_path", None) or "psql"
        self.database_user = kwds.get("postgres_database_user", None)
        self.database_host = kwds.get("postgres_database_host", None)
        self.database_port = kwds.get("postgres_database_port", None)
        self._kwds = kwds

    def sqlalchemy_url(self, identifier):
        """Return URL or form postgresql://username:password@localhost/mydatabase."""
        hostname = self.database_host or "localhost"
        if self.database_port:
            hostname += ":%s" % self.database_port
        return f"postgresql://{self.database_user}@{hostname}/{identifier}"

    def _psql_command_builder(self, *args):
        command_builder = _CommandBuilder(self.psql_path)
        # Print only tuples so output is easier to parse
        command_builder.append_command("--tuples-only")

        # Specify connection information
        if self.database_user:
            command_builder.append_command("--username", self.database_user)
        if self.database_host:
            command_builder.append_command("--host", self.database_host)
        if self.database_port:
            command_builder.append_command("--port", self.database_port)
        command_builder.append_command("-P", "pager=off")
        command_builder.extend_command(args)
        return command_builder


class _CommandBuilder:
    def __init__(self, *args):
        self.command = list(args)

    def append_command(self, *args_or_none):
        args_or_none = args_or_none or []
        for arg_or_none in args_or_none:
            if arg_or_none is not None:
                self.command.append(arg_or_none)

    def extend_command(self, args):
        for arg in args or []:
            self.append_command(arg)


__all__ = ("LocalPostgresDatabaseSource",)
