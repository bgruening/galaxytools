"""This modules describes the abstraction of a Galaxy profile.

This is a workspace with a specific default configuration and shed
tool setup. It is meant to be used with various serve commands.
"""

import json
import os
import shutil

import click
from galaxy.util.commands import which
from gxjobconfinit import (
    build_job_config,
    ConfigArgs,
)

from planemo.database import create_database_source
from planemo.galaxy.api import test_credentials_valid
from .config import DATABASE_LOCATION_TEMPLATE

PROFILE_OPTIONS_JSON_NAME = "planemo_profile_options.json"
ALREADY_EXISTS_EXCEPTION = "Cannot create profile with name [%s], directory [%s] already exists."


def profile_exists(ctx, profile_name, **kwds):
    """Return a truthy value iff the specified profile already exists."""
    profile_directory = _profile_directory(ctx, profile_name)
    return os.path.exists(profile_directory)


def list_profiles(ctx, **kwds):
    """Return a list of current profile names."""
    return os.listdir(ctx.galaxy_profiles_directory)


def delete_profile(ctx, profile_name, **kwds):
    """Delete profile with the specified name."""
    profile_directory = _profile_directory(ctx, profile_name)
    profile_options = _read_profile_options(profile_directory)
    profile_options, profile_options_path = _load_profile_to_json(ctx, profile_name)
    if profile_options["engine"] != "external_galaxy":
        database_type = profile_options.get("database_type")
        kwds["database_type"] = database_type
        if database_type != "sqlite":
            database_source = create_database_source(**kwds)
            database_identifier = _profile_to_database_identifier(profile_name)
            database_source.delete_database(
                database_identifier,
            )
    shutil.rmtree(profile_directory)


def create_profile(ctx, profile_name, **kwds):
    """Create a profile with the specified name."""
    engine_type = kwds.get("engine", "galaxy")
    profile_directory = _profile_directory(ctx, profile_name)
    if profile_exists(ctx, profile_name, **kwds):
        message = ALREADY_EXISTS_EXCEPTION % (profile_name, profile_directory)
        raise click.ClickException(message)

    os.makedirs(profile_directory)

    if engine_type == "docker_galaxy":
        create_for_engine = _create_profile_docker
    elif engine_type == "external_galaxy" or kwds.get("galaxy_url"):
        create_for_engine = _create_profile_external
    else:
        create_for_engine = _create_profile_local

    stored_profile_options = create_for_engine(ctx, profile_directory, profile_name, kwds)

    profile_options_path = _stored_profile_options_path(profile_directory)
    with open(profile_options_path, "w") as f:
        json.dump(stored_profile_options, f)


def _create_profile_docker(ctx, profile_directory, profile_name, kwds):
    export_directory = os.path.join(profile_directory, "export")
    os.makedirs(export_directory)
    return {
        "engine": "docker_galaxy",
    }


def _create_profile_local(ctx, profile_directory, profile_name, kwds):
    database_type = kwds.get("database_type", "auto")
    allow_sqlite_fallback = database_type == "auto"
    if database_type == "auto":
        if which("psql"):
            database_type = "postgres"
        elif which("docker"):
            database_type = "postgres_docker"
        elif which("singularity"):
            database_type = "postgres_singularity"
        else:
            database_type = "sqlite"

    if database_type not in ["sqlite", "postgres_singularity"]:
        database_source = create_database_source(**kwds)
        database_identifier = _profile_to_database_identifier(profile_name)
        try:
            database_source.create_database(
                database_identifier,
            )
        except RuntimeError:
            if allow_sqlite_fallback:
                # If postgres database creation fails (e.g., role doesn't exist, connection issues),
                # fall back to sqlite
                database_type = "sqlite"
            else:
                raise
        else:
            database_connection = database_source.sqlalchemy_url(database_identifier)
    elif database_type == "postgres_singularity":
        database_connection + database_source.sqlalchemy_url(database_identifier)
    if database_type == "sqlite":
        database_location = os.path.join(profile_directory, "galaxy.sqlite")
        database_connection = DATABASE_LOCATION_TEMPLATE % database_location

    return {
        "database_type": database_type,
        "database_connection": database_connection,
        "engine": "galaxy",
    }


def _create_profile_external(ctx, profile_directory, profile_name, kwds):
    url = kwds.get("galaxy_url")
    api_key = kwds.get("galaxy_admin_key") or kwds.get("galaxy_user_key")
    if test_credentials_valid(url=url, key=api_key, is_admin=kwds.get("galaxy_admin_key")):
        return {
            "galaxy_url": url,
            "galaxy_user_key": kwds.get("galaxy_user_key"),
            "galaxy_admin_key": kwds.get("galaxy_admin_key"),
            "engine": "external_galaxy",
        }
    else:
        raise ConnectionError("The credentials provided for an external Galaxy instance are not valid.")


def ensure_profile(ctx, profile_name, **kwds):
    """Ensure a Galaxy profile exists and return profile defaults."""
    if not profile_exists(ctx, profile_name, **kwds):
        available_profiles = list_profiles(ctx, **kwds)
        error_message = f"Profile '{profile_name}' does not exist."
        if available_profiles:
            error_message += f"\n\nAvailable profiles: {', '.join(available_profiles)}"
            error_message += f"\n\nTo create a new profile, use: planemo profile_create {profile_name}"
        else:
            error_message += (
                f"\n\nNo profiles found. To create a new profile, use: planemo profile_create {profile_name}"
            )
        raise click.UsageError(error_message)

    return _profile_options(ctx, profile_name, **kwds)


def create_alias(ctx, alias, obj, profile_name, **kwds):
    profile_options, profile_options_path = _load_profile_to_json(ctx, profile_name)

    if profile_options.get("aliases"):
        profile_options["aliases"][alias] = obj
    else:  # no aliases yet defined
        profile_options["aliases"] = {alias: obj}

    with open(profile_options_path, "w") as f:
        json.dump(profile_options, f)

    return 0


def list_alias(ctx, profile_name, **kwds):
    profile_options, _ = _load_profile_to_json(ctx, profile_name)
    return profile_options.get("aliases", {})


def delete_alias(ctx, alias, profile_name, **kwds):
    profile_options, profile_options_path = _load_profile_to_json(ctx, profile_name)
    if alias not in profile_options.get("aliases", {}):
        return 1
    else:
        del profile_options["aliases"][alias]

    with open(profile_options_path, "w") as f:
        json.dump(profile_options, f)

    return 0


def translate_alias(ctx, alias, profile_name):
    if not profile_name:
        return alias
    aliases = _load_profile_to_json(ctx, profile_name)[0].get("aliases", {})
    return aliases.get(alias, alias)


def _load_profile_to_json(ctx, profile_name):
    if not profile_exists(ctx, profile_name):
        raise click.ClickException("That profile does not exist. Create it with `planemo profile_create`")
    profile_directory = _profile_directory(ctx, profile_name)
    profile_options_path = _stored_profile_options_path(profile_directory)
    with open(profile_options_path) as f:
        profile_options = json.load(f)
    return profile_options, profile_options_path


def _profile_options(ctx, profile_name, **kwds):
    profile_directory = _profile_directory(ctx, profile_name)
    profile_options = _read_profile_options(profile_directory)

    if profile_options["engine"] == "docker_galaxy":
        engine_options = dict(export_directory=os.path.join(profile_directory, "export"))
    else:
        file_path = os.path.join(profile_directory, "files")
        shed_tool_path = os.path.join(profile_directory, "shed_tools")
        shed_tool_conf = os.path.join(profile_directory, "shed_tool_conf.xml")
        tool_dependency_dir = os.path.join(profile_directory, "deps")

        engine_options = dict(
            file_path=file_path,
            tool_dependency_dir=tool_dependency_dir,
            shed_tool_conf=shed_tool_conf,
            shed_tool_path=shed_tool_path,
            galaxy_brand=profile_name,
        )
    profile_options.update(engine_options)
    profile_options["galaxy_brand"] = profile_name
    return profile_options


def _profile_to_database_identifier(profile_name):
    char_lst = [c if c.isalnum() else "_" for c in profile_name]
    return "plnmoprof_%s" % "".join(char_lst)


def _read_profile_options(profile_directory):
    profile_options_path = _stored_profile_options_path(profile_directory)
    with open(profile_options_path) as f:
        profile_options = json.load(f)
    return profile_options


def _stored_profile_options_path(profile_directory):
    profile_options_path = os.path.join(profile_directory, PROFILE_OPTIONS_JSON_NAME)
    return profile_options_path


def _profile_directory(ctx, profile_name):
    return os.path.join(ctx.galaxy_profiles_directory, profile_name)


def initialize_job_config(ctx, profile_name, **kwds):
    profile_directory = _profile_directory(ctx, profile_name)
    job_config_path = os.path.join(profile_directory, "job_conf.yml")
    if os.path.exists(job_config_path):
        raise click.ClickException(f"File '{job_config_path}' already exists, exiting.")

    init_config = ConfigArgs.from_dict(**kwds)
    job_config = build_job_config(init_config)
    with open(job_config_path, "w") as f:
        f.write(job_config)

    config, profile_config_path = _load_profile_to_json(ctx, profile_name)
    config["job_config_file"] = os.path.abspath(profile_config_path)
    with open(profile_config_path, "w") as f:
        json.dump(config, f)
    return job_config_path


__all__ = (
    "create_profile",
    "delete_profile",
    "ensure_profile",
    "list_profiles",
)
