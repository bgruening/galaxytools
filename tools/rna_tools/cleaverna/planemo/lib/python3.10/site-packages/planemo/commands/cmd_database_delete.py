"""Module describing the planemo ``database_delete`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.database import create_database_source


@click.command("database_delete")
@options.database_identifier_argument()
@options.profile_database_options()
@options.docker_config_options()
@command_function
def cli(ctx, identifier, **kwds):
    """Delete a *development* database.

    Currently the only implementation is postgres which will be managed with
    ``psql``.

    Planemo ``database_`` commands make it very easy to create and destroy
    databases, therefore it should not be used for production data - and it
    should not even be connnected to a production database server. Planemo
    is intended for development purposes only.

    Planemo will assume that it can manage and access postgres databases
    without specifying a password. This can be accomplished by configuring
    postgres to not required a password for the planemo user or by specifying
    a password in a ``.pgpass`` file.

    Planemo can be configured to not require a password for the planemo user in
    the postgres configuration file ``pg_hba.conf`` (on Ubuntu/Debian linux
    distros this file is in /etc/postgresql/<postgres_version>/main/ directory).
    Adding the following lines to that file will allow planemo and Galaxy to
    access the databases without a password.

    \b
        # "local" is for Unix domain socket connections only
        local   all   all                    trust
        # IPv4 local connections:
        host    all   all    127.0.0.1/32    trust
        # IPv6 local connections:
        host    all   all    ::1/128         trust

    More information on the ``pg_hda.conf`` configuration file can be found at
    http://www.postgresql.org/docs/9.3/static/auth-pg-hba-conf.html.

    Information on ``.pgpass`` files can be found at at the following location:
    http://www.postgresql.org/docs/9.4/static/libpq-pgpass.html. In Ubuntu and
    Debian distros - a postgres user likely already exists and its password can
    be set by setting up a file ``~/.pgpass`` file with the following contents.

    \b
        *:*:*:postgres:<postgres_password>
    """
    create_database_source(**kwds).delete_database(identifier)
