"""General support infrastructure not tied to any particular test."""

import os
import random
import string
import unittest
from typing import (
    Any,
    Callable,
    Optional,
)

import requests

from bioblend.galaxy import GalaxyInstance

NO_GALAXY_MESSAGE = "Externally configured Galaxy required, but not found. Set BIOBLEND_GALAXY_URL and BIOBLEND_GALAXY_API_KEY to run this test."


def random_string(length: int = 8) -> str:
    return "".join(random.choice(string.ascii_lowercase) for _ in range(length))


def is_site_up(url: str) -> bool:
    try:
        response = requests.get(url, timeout=10)
        return response.status_code == 200
    except Exception:
        return False


def skip_unless_toolshed() -> Callable:
    """Decorate tests with this to skip the test if a URL for a ToolShed
    to run the tests is not provided.
    """
    if "BIOBLEND_TOOLSHED_URL" not in os.environ:
        return unittest.skip(
            "Externally configured ToolShed required, but not found. Set BIOBLEND_TOOLSHED_URL (e.g. to https://testtoolshed.g2.bx.psu.edu/ ) to run this test."
        )
    toolshed_url = os.environ["BIOBLEND_TOOLSHED_URL"]
    if not is_site_up(toolshed_url):
        return unittest.skip(f"Configured ToolShed [{toolshed_url}] appears to be down")
    return lambda f: f


def skip_unless_galaxy(min_release: Optional[str] = None) -> Callable:
    """Decorate tests with this to skip the test if Galaxy is not
    configured.
    """
    if min_release is not None:
        galaxy_release = os.environ.get("GALAXY_VERSION", None)
        if galaxy_release is not None and galaxy_release != "dev":
            if not galaxy_release.startswith("release_"):
                raise ValueError("The value of GALAXY_VERSION environment variable should start with 'release_'")
            if not min_release.startswith("release_"):
                raise Exception("min_release should start with 'release_'")
            if galaxy_release[8:] < min_release[8:]:
                return unittest.skip(f"Testing on Galaxy {galaxy_release}, but need {min_release} to run this test.")

    if "BIOBLEND_GALAXY_URL" not in os.environ:
        return unittest.skip(NO_GALAXY_MESSAGE)

    if "BIOBLEND_GALAXY_API_KEY" not in os.environ and "BIOBLEND_GALAXY_MASTER_API_KEY" in os.environ:
        galaxy_url = os.environ["BIOBLEND_GALAXY_URL"]
        galaxy_master_api_key = os.environ["BIOBLEND_GALAXY_MASTER_API_KEY"]
        gi = GalaxyInstance(galaxy_url, key=galaxy_master_api_key)

        if "BIOBLEND_GALAXY_USER_EMAIL" in os.environ:
            galaxy_user_email = os.environ["BIOBLEND_GALAXY_USER_EMAIL"]
        else:
            galaxy_user_email = f"{random_string()}@localhost.localdomain"

        galaxy_user_id = None
        for user in gi.users.get_users():
            if user["email"] == galaxy_user_email:
                galaxy_user_id = user["id"]
                break

        config = gi.config.get_config()
        if galaxy_user_id is None:
            if config.get("use_remote_user", False):
                new_user = gi.users.create_remote_user(galaxy_user_email)
            else:
                galaxy_user = galaxy_user_email.split("@", 1)[0]
                galaxy_password = random_string(20)

                # Create a new user
                new_user = gi.users.create_local_user(galaxy_user, galaxy_user_email, galaxy_password)
            galaxy_user_id = new_user["id"]

        if config["version_major"] >= "21.01":
            api_key = gi.users.get_or_create_user_apikey(galaxy_user_id)
        else:
            api_key = gi.users.get_user_apikey(galaxy_user_id)
            if not api_key or api_key == "Not available.":
                api_key = gi.users.create_user_apikey(galaxy_user_id)
        os.environ["BIOBLEND_GALAXY_API_KEY"] = api_key

    if "BIOBLEND_GALAXY_API_KEY" not in os.environ:
        return unittest.skip(NO_GALAXY_MESSAGE)

    return lambda f: f


def skip_unless_tool(tool_id: str) -> Callable:
    """Decorate a Galaxy test method as requiring a specific tool,
    skip the test case if the tool is unavailable.
    """

    def method_wrapper(method):
        def wrapped_method(has_gi, *args, **kwargs):
            tools = has_gi.gi.tools.get_tools()
            # In panels by default, so flatten out sections...
            tool_ids = [_["id"] for _ in tools]
            if tool_id not in tool_ids:
                raise unittest.SkipTest(f"Externally configured Galaxy instance requires tool {tool_id} to run test.")

            return method(has_gi, *args, **kwargs)

        # Must preserve method name so pytest can detect and report tests by
        # name.
        wrapped_method.__name__ = method.__name__
        return wrapped_method

    return method_wrapper


def get_abspath(path: str) -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), path))


def new_user_gi(admin_gi: GalaxyInstance) -> tuple[dict[str, Any], GalaxyInstance]:
    """Create a new local user and connect to the Galaxy instance as that user

    :param gi: GalaxyInstance object with admin privileges
    :return: a tuple with the new local user and the new GalaxyInstance object
      connected as that user
    """
    new_username = random_string()
    new_user_email = f"{new_username}@example.org"
    new_password = random_string(20)
    new_user = admin_gi.users.create_local_user(new_username, new_user_email, new_password)
    assert new_user["username"] == new_username
    assert new_user["email"] == new_user_email
    new_gi = GalaxyInstance(url=admin_gi.base_url, email=new_user_email, password=new_password)
    return new_user, new_gi
