"""Utilities for interacting with Github."""

import io
import os
import stat
import tarfile
import tempfile
import time
from distutils.dir_util import copy_tree
from pathlib import Path

import requests

from planemo import git
from planemo.io import (
    communicate,
    IS_OS_X,
)

try:
    import github

    has_github_lib = True
except ImportError:
    github = None
    has_github_lib = False

GH_VERSION = "1.5.0"

NO_GITHUB_DEP_ERROR = "Cannot use github functionality - PyGithub library not available."
FAILED_TO_DOWNLOAD_GH = "No gh executable available and it could not be installed."
DEFAULT_REMOTE_NAME = "planemo-remote"
SLEEP_BEFORE_RELEASE = int(os.environ.get("PLANEMO_SLEEP_BEFORE_RELEASE", 60))


def get_github_config(ctx, allow_anonymous=False):
    """Return a :class:`planemo.github_util.GithubConfig` for given configuration."""
    global_github_config = _get_raw_github_config(ctx)
    return GithubConfig(global_github_config, allow_anonymous=allow_anonymous)


def clone_fork_branch(ctx, target, path, remote_name=DEFAULT_REMOTE_NAME, **kwds):
    """Clone, fork, and branch a repository ahead of building a pull request."""
    git.checkout(ctx, target, path, branch=kwds.get("branch", None), remote="origin", from_branch="master")
    if kwds.get("fork"):
        try:
            fork(ctx, path, remote_name=remote_name, **kwds)
        except Exception:
            pass
    if "GITHUB_USER" in os.environ:
        # On CI systems fork doesn't add a local remote under circumstances I don't quite understand,
        # but that's probably linked to https://github.com/cli/cli/issues/2722
        cmd = [
            "git",
            "remote",
            "add",
            remote_name,
            f"https://github.com/{os.environ['GITHUB_USER']}/{os.path.basename(target)}",
        ]
        try:
            communicate(cmd, cwd=path)
        except RuntimeError:
            # Can add the remote only once
            pass
    return remote_name


def fork(ctx, path, remote_name=DEFAULT_REMOTE_NAME, **kwds):
    """Fork the target repository using ``gh``."""
    gh_path = ensure_gh(ctx, **kwds)
    gh_env = get_gh_env(ctx, path, **kwds)
    cmd = [gh_path, "repo", "fork", "--remote=true", "--remote-name", remote_name]
    communicate(cmd, cwd=path, env=gh_env)
    return remote_name


def get_or_create_repository(ctx, owner, repo, dry_run=True, **kwds):
    """Clones or creates a repository and returns path on disk"""
    target = os.path.realpath(tempfile.mkdtemp())
    remote_repo = f"https://github.com/{owner}/{repo}"
    try:
        ctx.log(f"Cloning {remote_repo}")
        git.clone(ctx, src=remote_repo, dest=target)
    except Exception:
        ctx.log(f"Creating repository {remote_repo}")
        target = create_repository(ctx, owner=owner, repo=repo, dest=target, dry_run=dry_run)
    return target


def create_repository(ctx, owner, repo, dest, dry_run, **kwds):
    gh_path = ensure_gh(ctx, **kwds)
    gh_env = get_gh_env(ctx, dry_run=dry_run, **kwds)
    cmd = [gh_path, "repo", "create", "-y", "--public", f"{owner}/{repo}"]
    if dry_run:
        "Would run command '{}'".format(" ".join(cmd))
        git.init(ctx, dest)
        return dest
    communicate(cmd, env=gh_env, cwd=dest)
    return os.path.join(dest, repo)


def rm_dir_contents(directory, ignore_dirs=(".git")):
    directory = Path(directory)
    for item in directory.iterdir():
        if item.name not in ignore_dirs:
            if item.is_dir():
                rm_dir_contents(item)
            else:
                item.unlink()


def add_dir_contents_to_repo(
    ctx, from_dir, target_dir, target_repository_path, version, dry_run, branch="main", notes=""
):
    ctx.log(f"From {from_dir} to {target_repository_path}")
    rm_dir_contents(target_repository_path)
    copy_tree(from_dir, target_repository_path)
    git.add(ctx, target_repository_path, target_repository_path)
    message = f"Update for version {version}"
    if notes:
        message += f"\n{notes}"
    git.commit(ctx, repo_path=target_repository_path, message=message)
    if not dry_run:
        git.push(ctx, target_repository_path, to="origin", branch=branch)


def assert_new_version(ctx, version, owner, repo):
    remote_repo = f"https://github.com/{owner}/{repo}"
    try:
        tags_and_versions = git.ls_remote(ctx, remote_repo=remote_repo)
        if f"refs/tags/v{version}" in tags_and_versions or f"refs/tags/{version}" in tags_and_versions:
            raise Exception(f"Version '{version}' for {owner}/{repo} exists already. Please change the version.")
    except RuntimeError:
        # repo doesn't exist
        pass


def changelog_in_repo(target_repository_path):
    """
    Parse the changelog for the current version.

    The changelog should be written following the conventions in https://keepachangelog.com/en/0.3.0/.
    This means we will assume that individual versions are announced using `##`.
    """
    changelog = []
    for path in os.listdir(target_repository_path):
        if "changelog.md" in path.lower():
            version_header_seen = False
            version_header_prefix = "## "
            top_level_header_seen = False
            top_level_header_prefix = "# "
            with open(os.path.join(target_repository_path, path)) as changelog_fh:
                for line in changelog_fh:
                    if line.startswith(top_level_header_prefix):
                        top_level_header_seen = True
                        continue
                    if line.startswith(version_header_prefix):
                        if version_header_seen:
                            return "\n".join(changelog[:-1])
                        else:
                            version_header_seen = True
                    if top_level_header_seen and not line or not version_header_seen:
                        continue
                    changelog.append(line.rstrip())
    return "\n".join(changelog).rstrip()


def create_release(ctx, from_dir, target_dir, owner, repo, version, dry_run, branch, notes="", **kwds):
    assert_new_version(ctx, version, owner=owner, repo=repo)
    target_repository_path = get_or_create_repository(ctx, owner=owner, repo=repo, dry_run=dry_run)
    add_dir_contents_to_repo(
        ctx, from_dir, target_dir, target_repository_path, version=version, dry_run=dry_run, branch=branch, notes=notes
    )
    gh_path = ensure_gh(ctx, **kwds)
    gh_env = get_gh_env(ctx, dry_run=dry_run, **kwds)
    cmd = [
        gh_path,
        "release",
        "-R",
        f"{owner}/{repo}",
        "create",
        f"v{version}",
        "--title",
        str(version),
    ]
    cmd.extend(["--notes", notes or changelog_in_repo(target_repository_path)])
    if not dry_run:
        # For new repositories dockstore needs a bit of time to register the new workflow.
        time.sleep(SLEEP_BEFORE_RELEASE)
        communicate(cmd, env=gh_env)
    else:
        ctx.log("Would run command '{}'".format(" ".join(cmd)))


def pull_request(ctx, path, message=None, repo=None, **kwds):
    """Create a pull request against the origin of the path using ``gh``."""
    gh_path = ensure_gh(ctx, **kwds)
    gh_env = get_gh_env(ctx, path, **kwds)
    cmd = [gh_path, "pr", "create"]
    if message is None:
        cmd.append("--fill")
    else:
        lines = message.splitlines()
        cmd.extend(["--title", lines[0]])
        if len(lines) > 1:
            cmd.extend(["--body", "\n".join(lines[1:])])
    if repo:
        cmd.extend(["--repo", repo])
    communicate(cmd, env=gh_env)


def get_gh_env(ctx, path=None, dry_run=False, **kwds):
    """Return a environment dictionary to run gh with given user and repository target."""
    if path is None:
        env = {}
    else:
        env = git.git_env_for(path).copy()
    if not dry_run:
        github_config = _get_raw_github_config(ctx)
        if github_config is not None:
            if "access_token" in github_config:
                env["GITHUB_TOKEN"] = github_config["access_token"]

    return env


def ensure_gh(ctx, **kwds):
    """Ensure gh is available for planemo

    This method will ensure ``gh`` is installed at the correct version.

    For more information on ``gh`` checkout https://cli.github.com/
    """
    planemo_gh_path = os.path.join(ctx.workspace, f"gh-{GH_VERSION}")
    if not os.path.exists(planemo_gh_path):
        _try_download_gh(planemo_gh_path)

    if not os.path.exists(planemo_gh_path):
        raise Exception(FAILED_TO_DOWNLOAD_GH)

    return planemo_gh_path


def _try_download_gh(planemo_gh_path):
    link = _gh_link()
    path = Path(planemo_gh_path)
    resp = requests.get(link)
    with tarfile.open(fileobj=io.BytesIO(resp.content)) as tf, path.open("wb") as outfile:
        for member in tf.getmembers():
            if member.name.endswith("bin/gh"):
                outfile.write(tf.extractfile(member).read())
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _get_raw_github_config(ctx):
    """Return a :class:`planemo.github_util.GithubConfig` for given configuration."""
    if "github" not in ctx.global_config:
        if "GITHUB_TOKEN" in os.environ:
            return {
                "access_token": os.environ["GITHUB_TOKEN"],
            }
    if "github" not in ctx.global_config:
        raise Exception("github account not found in planemo config and GITHUB_TOKEN environment variables unset")
    return ctx.global_config["github"]


class GithubConfig:
    """Abstraction around a Github account.

    Required to use ``github`` module methods that require authorization.
    """

    def __init__(self, config, allow_anonymous=False):
        if not has_github_lib:
            raise Exception(NO_GITHUB_DEP_ERROR)
        if "access_token" not in config:
            if not allow_anonymous:
                raise Exception("github authentication unavailable")
            github_object = github.Github()
        else:
            github_object = github.Github(config["access_token"])
        self._github = github_object


def _gh_link():
    if IS_OS_X:
        template_link = "https://github.com/cli/cli/releases/download/v%s/gh_%s_macOS_amd64.tar.gz"
    else:
        template_link = "https://github.com/cli/cli/releases/download/v%s/gh_%s_linux_amd64.tar.gz"
    return template_link % (GH_VERSION, GH_VERSION)


def publish_as_gist_file(ctx, path, name="index"):
    """Publish a gist.

    More information on gists at http://gist.github.com/.
    """
    github_config = get_github_config(ctx, allow_anonymous=False)
    user = github_config._github.get_user()
    with open(path) as fh:
        content = fh.read()
    content_file = github.InputFileContent(content)
    gist = user.create_gist(False, {name: content_file})
    return gist.files[name].raw_url


def get_repository_object(ctx, name):
    github_object = get_github_config(ctx, allow_anonymous=True)
    return github_object._github.get_repo(name)


__all__ = (
    "add_dir_contents_to_repo",
    "clone_fork_branch",
    "create_release",
    "ensure_gh",
    "fork",
    "get_github_config",
    "get_gh_env",
    "get_or_create_repository",
    "publish_as_gist_file",
)
