"""Module describes a :class:`DatabaseSource` for managed, dockerized postgres databases."""

import os
import subprocess
import time
from tempfile import mkdtemp

from galaxy.util.commands import execute

from planemo.io import info
from .interface import DatabaseSource
from .postgres import ExecutesPostgresSqlMixin

DEFAULT_CONTAINER_NAME = "planemopostgres"
DEFAULT_POSTGRES_DATABASE_NAME = "galaxy"
DEFAULT_POSTGRES_USER = "galaxy"
DEFAULT_POSTGRES_PASSWORD = "mysecretpassword"
DEFAULT_POSTGRES_PORT_EXPOSE = 5432
DEFAULT_DOCKERIMAGE = "postgres:14.2-alpine3.15"
DEFAULT_SIF_NAME = "postgres_14_2-alpine3_15.sif"

DEFAULT_CONNECTION_STRING = f"postgresql://{DEFAULT_POSTGRES_USER}:{DEFAULT_POSTGRES_PASSWORD}@localhost:{DEFAULT_POSTGRES_PORT_EXPOSE}/{DEFAULT_POSTGRES_DATABASE_NAME}"


def start_postgres_singularity(
    singularity_path,
    container_instance_name,
    database_location,
    databasename=DEFAULT_POSTGRES_DATABASE_NAME,
    user=DEFAULT_POSTGRES_USER,
    password=DEFAULT_POSTGRES_PASSWORD,
    **kwds,
):
    info(f"Postgres database stored at: {database_location}")
    pgdata_path = os.path.join(database_location, "pgdata")
    pgrun_path = os.path.join(database_location, "pgrun")

    if not os.path.exists(pgdata_path):
        os.makedirs(pgdata_path)
    if not os.path.exists(pgrun_path):
        os.makedirs(pgrun_path)

    version_file = os.path.join(pgdata_path, "PG_VERSION")
    if not os.path.exists(version_file):
        # Run container for a short while to initialize the database
        # The database will not be initilizaed during a
        # "singularity instance start" command
        init_database_command = [
            singularity_path,
            "run",
            "-B",
            f"{pgdata_path}:/var/lib/postgresql/data",
            "-B",
            f"{pgrun_path}:/var/run/postgresql",
            "-e",
            "-C",
            "--env",
            f"POSTGRES_DB={databasename}",
            "--env",
            f"POSTGRES_USER={user}",
            "--env",
            f"POSTGRES_PASSWORD={password}",
            "--env",
            "POSTGRES_INITDB_ARGS='--encoding=UTF-8'",
            f"docker://{DEFAULT_DOCKERIMAGE}",
        ]
        info(f"Initilizing postgres database in folder: {pgdata_path}")
        process = subprocess.Popen(init_database_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Give the container time to initialize the database
        for _ in range(10):
            if os.path.exists(version_file):
                break
            time.sleep(5)
            info("Waiting for the postgres database to initialize.")
        else:
            raise Exception("Failed to initialize the postgres database.")
        time.sleep(10)
        process.terminate()

    # Start the singularity instance, assumes the database is
    # already initialized since the entrypoint will not be run
    # when starting a instance of the container.
    run_command = [
        singularity_path,
        "instance",
        "start",
        "-B",
        f"{pgdata_path}:/var/lib/postgresql/data",
        "-B",
        f"{pgrun_path}:/var/run/postgresql",
        "-e",
        "-C",
        f"docker://{DEFAULT_DOCKERIMAGE}",
        container_instance_name,
    ]
    info(f"Starting singularity instance named: {container_instance_name}")
    execute(run_command)


def stop_postgress_singularity(container_instance_name, **kwds):
    info(f"Stopping singularity instance named: {container_instance_name}")
    execute(["singularity", "instance", "stop", container_instance_name])


class SingularityPostgresDatabaseSource(ExecutesPostgresSqlMixin, DatabaseSource):
    """
    Postgres database running inside a Singularity container. Should be used with
    "with" statements to automatically start and stop the container.
    """

    def __init__(self, **kwds):
        """Construct a postgres database source from planemo configuration."""

        self.singularity_path = "singularity"
        self.database_user = DEFAULT_POSTGRES_USER
        self.database_password = DEFAULT_POSTGRES_PASSWORD
        self.database_host = "localhost"  # TODO: Make docker host
        self.database_port = DEFAULT_POSTGRES_PORT_EXPOSE
        if "postgres_storage_location" in kwds and kwds["postgres_storage_location"] is not None:
            self.database_location = kwds["postgres_storage_location"]
        else:
            self.database_location = os.path.join(mkdtemp(suffix="_planemo_postgres_db"))
        self.container_instance_name = f"{DEFAULT_CONTAINER_NAME}-{int(time.time() * 1000000)}"
        self._kwds = kwds

    def __enter__(self):
        start_postgres_singularity(
            singularity_path=self.singularity_path,
            database_location=self.database_location,
            user=self.database_user,
            password=self.database_password,
            container_instance_name=self.container_instance_name,
            **self._kwds,
        )

    def __exit__(self, exc_type, exc_value, traceback):
        stop_postgress_singularity(self.container_instance_name)

    def sqlalchemy_url(self, identifier):
        """Return URL or form postgresql://username:password@localhost/mydatabase."""
        return "postgresql://%s:%s@%s:%d/%s" % (
            self.database_user,
            self.database_password,
            self.database_host,
            self.database_port,
            identifier,
        )


__all__ = ("SingularityPostgresDatabaseSource",)
