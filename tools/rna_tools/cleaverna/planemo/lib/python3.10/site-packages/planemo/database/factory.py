"""Create a DatabaseSource from supplied planemo configuration."""

from galaxy.util.commands import which

from .interface import DatabaseSource
from .postgres import LocalPostgresDatabaseSource
from .postgres_docker import DockerPostgresDatabaseSource
from .postgres_singularity import SingularityPostgresDatabaseSource


def create_database_source(**kwds) -> DatabaseSource:
    """Return a :class:`planemo.database.interface.DatabaseSource` for configuration."""
    database_type = kwds.get("database_type", "auto")
    if database_type == "auto":
        if which("psql"):
            database_type = "postgres"
        elif which("docker"):
            database_type = "postgres_docker"
        elif which("singularity"):
            database_type = "postgres_singularity"
        else:
            raise Exception("Cannot find executables for psql or docker, cannot configure a database source.")

    if database_type == "postgres":
        return LocalPostgresDatabaseSource(**kwds)
    elif database_type == "postgres_docker":
        return DockerPostgresDatabaseSource(**kwds)
    elif database_type == "postgres_singularity":
        return SingularityPostgresDatabaseSource(**kwds)
    # TODO
    # from .sqlite import SqliteDatabaseSource
    # elif database_type == "sqlite":
    #     return SqliteDatabaseSource(**kwds)
    else:
        raise Exception("Unknown database type [%s]." % database_type)


__all__ = ("create_database_source",)
