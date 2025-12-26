"""Interface over bioblend and direct access to ToolShed API via requests."""

import json
from typing import Optional

from bioblend.toolshed import ToolShedInstance
from galaxy.util import unicodify

from planemo.io import untar_to

REPOSITORY_DOWNLOAD_TEMPLATE = "%s/repository/download?repository_id=%s&changeset_revision=default&file_type=gz"


def tool_shed_instance(
    url: str, key: Optional[str] = None, email: Optional[str] = None, password: Optional[str] = None
) -> ToolShedInstance:
    return ToolShedInstance(url=url, key=key, email=email, password=password)


def find_repository(tsi, owner, name):
    """Find repository information for given owner and repository
    name.
    """
    repos = tsi.repositories.get_repositories(owner=owner, name=name)

    def matches(r):
        return r["owner"] == owner and r["name"] == name

    matching_repos = list(filter(matches, repos))
    if not matching_repos:
        return None
    else:
        return matching_repos[0]


def latest_installable_revision(tsi, repository_id):
    info = tsi.repositories.show_repository(repository_id)
    owner = info["owner"]
    name = info["name"]
    revisions = tsi.repositories.get_ordered_installable_revisions(name, owner)
    if len(revisions) == 0:
        msg = f"Failed to find installable revisions for [{owner}, {name}]."
        raise Exception(msg)
    else:
        return revisions[-1]


def username(tsi):
    """Fetch current username from shed given API key/auth."""
    user = _user(tsi)
    return user["username"]


def api_exception_to_message(e):
    """Convert API exception to human digestable error message - parsing
    out the shed generate message if possible.
    """
    message = unicodify(e)
    if hasattr(e, "read"):
        message = e.read()
        try:
            # Galaxy passes nice JSON messages as their errors, which bioblend
            # blindly returns. Attempt to parse those.
            upstream_error = json.loads(message)
            message = upstream_error["err_msg"]
        except Exception:
            pass
    return message


def find_category_ids(tsi, categories):
    """Translate human readable category names into their associated IDs."""
    category_list = tsi.categories.get_categories()

    category_ids = []
    for cat in categories:
        matching_cats = [x for x in category_list if x["name"] == cat]
        if not matching_cats:
            message = "Failed to find category %s" % cat
            raise Exception(message)
        category_ids.append(matching_cats[0]["id"])
    return category_ids


def download_tar(tsi, repo_id, destination, to_directory):
    base_url = tsi.base_url.rstrip("/")
    download_url = REPOSITORY_DOWNLOAD_TEMPLATE % (base_url, repo_id)
    if to_directory:
        tar_args = ["-xzf", "-", "--strip-components=1"]
        untar_to(download_url, tar_args=tar_args, dest_dir=destination)
    else:
        untar_to(download_url, path=destination)


def _user(tsi):
    """Fetch user information from the ToolShed API for given
    key.
    """
    # TODO: this should be done with an actual bioblend method,
    # see https://github.com/galaxyproject/bioblend/issues/130.
    response = tsi.make_get_request(tsi.url + "/users")
    return response.json()[0]
