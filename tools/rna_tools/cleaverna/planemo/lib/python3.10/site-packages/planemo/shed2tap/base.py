import os
import subprocess
import sys
import tarfile
import zipfile
from ftplib import all_errors as FTPErrors  # tuple of exceptions
from typing import List
from urllib.error import URLError
from urllib.request import urlretrieve
from xml.etree import ElementTree

from galaxy.util import unicodify

TOOLSHED_MAP = {
    "toolshed": "https://toolshed.g2.bx.psu.edu",
    "testtoolshed": "https://testtoolshed.g2.bx.psu.edu",
}


class Dependencies:
    """Base class for parsing Tool Shed dependency files."""

    def __init__(
        self,
        dependencies_file,
        repo=None,
        package_factory=None,
    ):
        if package_factory is None:
            package_factory = BasePackage
        self.repo = repo
        self.root = ElementTree.parse(dependencies_file).getroot()
        packages = []
        dependencies = []
        package_els = self.root.findall("package")
        assert package_els is not None
        for package_el in package_els:
            install_els = package_el.findall("install")
            readme_els = package_el.findall("readme")
            if len(readme_els) > 0:
                readme = readme_els[0].text
            else:
                readme = None
            assert len(install_els) in (0, 1)
            if len(install_els) == 1:
                install_el = install_els[0]
                package = package_factory(self, package_el, install_el, readme=readme)
                packages.append(package)
            else:
                repository_el = package_el.find("repository")
                if repository_el is None:
                    message = f"no repository in package el for {repo}"
                    raise AssertionError(message)
                dependency = Dependency(self, package_el, repository_el)
                dependencies.append(dependency)

        self.packages = packages
        self.dependencies = dependencies

    def single_package(self):
        return len(self.packages) == 1

    def __repr__(self):
        return f"Dependencies[for_repo={self.repo}]"


class Repo:
    def __init__(self, **kwds):
        for key, value in kwds.items():
            setattr(self, key, value)

    def recipe_base_name(self):
        owner = self.owner.replace("-", "")
        name = self.name
        name = name.replace("_", "").replace("-", "")
        base = f"{owner}_{name}"
        return base

    @staticmethod
    def from_xml(elem):
        tool_shed_url = elem.attrib.get("toolshed", None)
        if tool_shed_url and ("testtoolshed" in tool_shed_url):
            prefix = "testtoolshed"
        else:
            prefix = "toolshed"
        prior = elem.attrib.get("prior_installation_required", False)
        return Repo(
            prefix=prefix,
            name=elem.attrib["name"],
            owner=elem.attrib["owner"],
            tool_shed_url=tool_shed_url,
            changeset_revision=elem.attrib.get("changeset_revision", None),
            prior_installation_required=prior,
        )

    @staticmethod
    def from_api(prefix, repo_json):
        return Repo(
            prefix=prefix,
            name=repo_json["name"],
            owner=repo_json["owner"],
            tool_shed_url=TOOLSHED_MAP[prefix],
        )

    def get_file(self, path):
        try:
            url = f"{self.tool_shed_url}/repos/{self.owner}/{self.name}/raw-file/tip/{path}"
            path, headers = urlretrieve(url)
            return path
        except Exception as e:
            print(e)
            return None

    def __repr__(self):
        return f"Repository[name={self.name},owner={self.owner}]"


class Dependency:
    def __init__(self, dependencies, package_el, repository_el):
        self.dependencies = dependencies
        self.package_el = package_el
        self.repository_el = repository_el
        self.repo = Repo.from_xml(repository_el)

    def __repr__(self):
        return (
            f"Dependency[package_name={self.package_el.attrib['name']},version={self.package_el.attrib['version']},"
            f"dependent_package={self.repository_el.attrib['name']}]"
        )


class BasePackage:
    def __init__(self, dependencies, package_el, install_el, readme):
        self.dependencies = dependencies
        self.package_el = package_el
        self.install_el = install_el
        self.readme = readme
        self.all_actions = self.get_all_actions()
        self.no_arch_option = self.has_no_achitecture_install()

    def get_all_actions(self):
        action_or_group = self.install_el[0]
        parsed_actions = []
        if action_or_group.tag == "actions":
            parsed_actions.append(self.parse_actions(action_or_group))
        elif action_or_group.tag == "actions_group":
            actions_els = action_or_group.findall("actions")
            assert actions_els is not None
            for actions in actions_els:
                parsed_actions.append(self.parse_actions(actions))
            action_els = action_or_group.findall("action")
            assert action_els is not None
            for action in action_els:
                for parsed_a in parsed_actions:
                    parsed_a.actions.append(self.parse_action(action))
        return parsed_actions

    def has_no_achitecture_install(self):
        all_actions = self.all_actions
        if len(all_actions) < 2:
            return False
        else:
            last_action = all_actions[-1]
            return (not last_action.architecture) and (not last_action.os)

    def has_explicit_set_environments(self):
        all_actions = self.all_actions
        for actions in all_actions:
            for action in actions.actions:
                if action.explicit_variables:
                    return True
        return False

    def has_multiple_set_environments(self):
        all_actions = self.all_actions
        for actions in all_actions:
            count = 0
            for action in actions.actions:
                if action.explicit_variables:
                    count += 1
            if count > 1:
                return True
        return False

    def parse_actions(self, actions):
        os = actions.attrib.get("os", None)
        architecture = actions.get("architecture", None)
        action_els = actions.findall("action")
        assert action_els is not None
        parsed_actions = list(map(self.parse_action, action_els))
        action_packages = []
        for package in actions.findall("package"):
            action_packages.append(self.parse_action_package(package))
        return Actions(parsed_actions, os, architecture, action_packages)

    def parse_action_package(self, elem):
        name = elem.attrib["name"]
        version = elem.attrib["version"]
        repo = Repo.from_xml(elem.find("repository"))
        return ActionPackage(name, version, repo)

    def parse_action(self, action):
        return BaseAction.from_elem(action, package=self)

    def __repr__(self):
        actions = self.all_actions
        return (
            f"Install[name={self.package_el.attrib['name']},version={self.package_el.attrib['version']},"
            f"dependencies={self.dependencies},actions={actions}]"
        )


class Actions:
    def __init__(self, actions, os=None, architecture=None, action_packages=[]):
        self.os = os
        self.architecture = architecture
        self.actions = actions or []
        self.action_packages = action_packages

    def first_download(self):
        for action in self.actions:
            if action.action_type in ["download_by_url", "download_file"]:
                return action
        return None

    def downloads(self):
        actions = []
        for action in self.actions:
            if action.action_type in ["download_by_url", "download_file"]:
                actions.append(action)
        return actions

    def __repr__(self):
        platform = ""
        if self.os or self.architecture:
            platform = f"os={self.os},arch={self.architecture},"
        return f"Actions[{platform}{map(str, self.actions)}]"

    def _indent_extend(self, target, new_entries, indent="    "):
        for line in new_entries:
            target.append(indent + line)

    def to_bash(self):
        # Use self.os.title() to match "Linux" or "Darwin" in bash where case matters:
        if self.os and self.architecture:
            condition = f'("{self.os.title()}" == `uname`) && ("{self.architecture}" == `arch`)'
        elif self.os:
            condition = f'"{self.os.title()}" == `uname`'
        elif self.architecture:
            condition = f'"{self.architecture}" == `arch`'
        else:
            condition = None

        install_cmds = []
        env_cmds = []

        if condition:
            # Conditional actions block
            install_cmds = [
                "#" + "-" * 60,
                f"if [[ $specifc_action_done == 0 && {condition} ]]",
                "then",
                f'    echo "Platform-specific action for os={self.os}, arch={self.architecture}"',
            ]
            env_cmds = install_cmds[:]
            # TODO - Refactor block indentation?
            for action in self.actions:
                i_cmds, e_cmds = action.to_bash()
                self._indent_extend(install_cmds, i_cmds)
                self._indent_extend(env_cmds, e_cmds)
            # If we run the action, do not want to run any later actions!
            install_cmds.extend(["    specifc_action_done=1", "fi"])
            env_cmds.extend(["    specifc_action_done=1", "fi"])
        else:
            # Non-specific default action...
            install_cmds = [
                "#" + "-" * 60,
                "if [[ $specifc_action_done == 0 ]]",
                "then",
                '    echo "Non-platform-specific actions"',
            ]
            env_cmds = install_cmds[:]
            for action in self.actions:
                i_cmds, e_cmds = action.to_bash()
                self._indent_extend(install_cmds, i_cmds)
                self._indent_extend(env_cmds, e_cmds)
            install_cmds.append("fi")
            env_cmds.append("fi")
        return install_cmds, env_cmds


class ActionPackage:
    def __init__(self, name, version, repo):
        self.name = name
        self.version = version
        self.repo = repo


class BaseAction:
    _keys: List[str] = []
    action_type: str

    def __repr__(self):
        return f"Action[type={self.action_type}]"

    def same_as(self, other):
        if self._keys != other._keys:
            return False
        else:
            for key in self._keys:
                if getattr(self, key) != getattr(other, key):
                    return False

            return True

    def parse_action_repo(self, elem):
        repo_elem = elem.find("repository")
        repo = Repo.from_xml(repo_elem)
        self.repo = repo

    def parse_package_elems(self, elem):
        package_els = elem.findall("package")
        packages = []
        assert package_els is not None
        for package_el in package_els:
            packages.append(package_el.text)
        self.packages = packages

    @classmethod
    def from_elem(cls, elem, package):
        type = elem.attrib["type"]
        action_class = actions_by_type[type]
        return action_class(elem)

    def to_bash(self):
        """Return lists of bash shell commands to execute this action.

        This method is be implemented by each sub-class, and will
        return two list of strings (for ``dep_install.sh`` and
        ``env.sh`` respectively).
        """
        raise NotImplementedError(f"No to_bash defined for {self!r}")


def _tar_folders(filename):
    with tarfile.open(filename, "r", errorlevel=0) as archive:
        folders = set()
        for i in archive.getmembers():
            if i.isdir():
                folders.add(i.name.rstrip("/"))
            else:
                folders.add(os.path.split(i.name)[0])
        return list(folders)


def _zip_folders(filename):
    archive = zipfile.ZipFile(filename, "r")
    return list({i.filename.rstrip("/") for i in archive.infolist() if i.filename.endswith("/")})


def _common_prefix(folders):
    common_prefix = ""
    if len(folders) == 1:
        common_prefix = list(folders)[0]
    else:
        common_prefix = os.path.commonprefix(folders)
        assert not os.path.isabs(common_prefix), folders
    return common_prefix


def _cache_download(url, filename, sha256sum=None):
    """Returns local path to cached copy of URL using given filename."""
    cache = os.environ.get("DOWNLOAD_CACHE", "./download_cache/")
    # TODO - expose this as a command line option

    if not os.path.isdir(cache):
        os.mkdir(cache)

    local = os.path.join(cache, filename)

    if not os.path.isfile(local):
        # Must download it...
        try:
            # TODO - log this nicely...
            sys.stderr.write(f"Downloading {url} to {local!r}\n")
            urlretrieve(url, local)
        except URLError:
            # Most likely server is down, could be bad URL in XML action:
            raise RuntimeError(f"Unable to download {url}")
        except FTPErrors:
            # Most likely server is down, could be bad URL in XML action:
            raise RuntimeError(f"Unable to download {url}")

        # Verifying the checksum is slow, only do this on a fresh
        # download. Assume locally cached files are already OK.
        if sha256sum:
            # TODO - log this nicely...
            sys.stderr.write(f"Verifying checksum for {filename}\n")
            filehash = subprocess.check_output(["shasum", "-a", "256", local])[0:64].strip()
            filehash = unicodify(filehash)
            if filehash != sha256sum:
                raise RuntimeError(f"Checksum failure for {local}, got {filehash!r} but wanted {sha256sum!r}")

    return local


def _determine_compressed_file_folder(url, downloaded_filename, target_filename=None, sha256sum=None):
    """Determine how to decompress the file & its directory structure.

    Returns a list of shell commands. Consider this example where the
    folder to change to cannot be guessed from the tar-ball filename:

        $ curl -o "ncbi-blast-2.2.30+-ia32-linux.tar.gz" \
        "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-ia32-linux.tar.gz"
        $ tar -zxvf ncbi-blast-2.2.30+-ia32-linux.tar.gz
        $ cd ncbi-blast-2.2.30+

    Here it would return:

        ['tar -zxvf ncbi-blast-2.2.30+-ia32-linux.tar.gz', 'cd ncbi-blast-2.2.30+']

    If not cached, this function will download the file to the
    $DOWNLOAD_CACHE folder, and then open it / decompress it in
    order to find common folder prefix used.  This will also verify
    how to decompress the file, and the checksum if given.
    """
    answer = []

    local = _cache_download(url, downloaded_filename, sha256sum)

    if not target_filename:
        target_filename = downloaded_filename

    if tarfile.is_tarfile(local):
        folders = _tar_folders(local)
        if target_filename.endswith((".tar.gz", ".tgz")):
            answer.append(f"tar -zxvf {target_filename}")
        elif target_filename.endswith(".tar.bz2"):
            answer.append(f"tar -jxvf {target_filename}")
        elif target_filename.endswith(".tar"):
            answer.extend(f"tar -xvf {target_filename}")
        else:
            # Quite possibly this file doesn't need decompressing,
            # but until we've tested lots of real world tool_dependencies.xml
            # files I'd like to check these cases to confirm this.
            raise NotImplementedError(f"How to decompress tar file {target_filename}?")
    elif zipfile.is_zipfile(local):
        if target_filename.endswith(".jar"):
            # Do not decompress!
            return answer
        folders = _zip_folders(local)
        answer.append(f"unzip {target_filename}")
    elif target_filename.endswith(".dmg"):
        # Do not decompress!
        return answer
    else:
        # No compression? Leave as it is?
        raise NotImplementedError(f"What kind of compression is {local} using?")

    common_prefix = _common_prefix(folders)
    if common_prefix:
        answer.append(f'cd "{common_prefix}"')

    return answer


def _commands_and_downloaded_file(url, target_filename=None, sha256sum=None):
    # We preserve the filename from the URL in the cache.
    # i.e. We do NOT use the target_filename in the cache.
    # This because some Galaxy recipes normalise platform specific downloads
    # to use a single target filename, which would therefore break checksums etc
    # e.g. tests/data/repos/package_1/tool_dependencies.xml
    downloaded_filename = os.path.split(url)[-1]
    if "?" in downloaded_filename:
        downloaded_filename = downloaded_filename[: downloaded_filename.index("?")]
    if "#" in downloaded_filename:
        downloaded_filename = downloaded_filename[: downloaded_filename.index("#")]

    if not target_filename:
        target_filename = downloaded_filename

    # Curl is present on Mac OS X, can we assume it will be on Linux?
    # Cannot assume that wget will be on Mac OS X.
    answer = [
        f'if [[ -f "{target_filename}" ]]',
        "then",
        f'    echo "Reusing existing {target_filename}"',
        f'elif [[ -f "$DOWNLOAD_CACHE/{downloaded_filename}" ]]',
        "then",
        f'    echo "Reusing cached {downloaded_filename}"',
        f'    cp "$DOWNLOAD_CACHE/{downloaded_filename}" "{target_filename}"',
        "else",
        f'    echo "Downloading {downloaded_filename}"',
        f'    curl -L -o "$DOWNLOAD_CACHE/{downloaded_filename}" "{url}"',
        f'    cp "$DOWNLOAD_CACHE/{downloaded_filename}" "{target_filename}"',
    ]
    if sha256sum:
        # This is inserted into the if-else for a fresh download only.
        # Note double space between checksum and filename:
        answer.append(f'    echo "{sha256sum}  {target_filename}" | shasum -a 256 -c -')
    answer.append("fi")

    return answer, downloaded_filename


def _commands_to_download_and_extract(url, target_filename=None, sha256sum=None):
    answer, downloaded_filename = _commands_and_downloaded_file(url, target_filename, sha256sum)
    # Now should we unpack the tar-ball etc?
    answer.extend(_determine_compressed_file_folder(url, downloaded_filename, target_filename, sha256sum))
    return answer, []


class DownloadByUrlAction(BaseAction):
    action_type = "download_by_url"
    _keys = ["url"]

    def __init__(self, elem):
        self.url = elem.text.strip()
        assert self.url
        self.sha256sum = elem.attrib.get("sha256sum", None)
        self.target_filename = elem.attrib.get("target_filename", None)

    def to_bash(self):
        # See class DownloadByUrl in Galaxy,
        # lib/tool_shed/galaxy_install/tool_dependencies/recipe/step_handler.py
        return _commands_to_download_and_extract(self.url, self.target_filename, self.sha256sum)


class DownloadFileAction(BaseAction):
    action_type = "download_file"
    _keys = ["url", "extract"]

    def __init__(self, elem):
        self.url = elem.text.strip()
        self.extract = asbool(elem.attrib.get("extract", False))
        self.sha256sum = elem.attrib.get("sha256sum", None)
        self.target_filename = elem.attrib.get("target_filename", None)

    def to_bash(self):
        if self.extract:
            return _commands_to_download_and_extract(self.url, self.target_filename, self.sha256sum)
        else:
            commands, downloaded_file = _commands_and_downloaded_file(self.url, self.target_filename, self.sha256sum)
            return commands, []


class DownloadBinary(BaseAction):
    action_type = "download_binary"
    _keys = ["url_template", "target_directory"]

    def __init__(self, elem):
        self.url_template = elem.text
        assert self.url_template
        self.target_directory = elem.get("target_directory", None)


class ShellCommandAction(BaseAction):
    action_type = "shell_command"
    _keys = ["command"]

    def __init__(self, elem):
        self.command = elem.text

    def to_bash(self):
        # Galaxy would run each action from the same temp
        # working directory - possible that tool_dependencies.xml
        # shell_command could change $PWD so reset this:
        return ["pushd . > /dev/null", self.command, "popd > /dev/null"], []


class TemplateShellCommandAction(BaseAction):
    action_type = "template_command"
    _keys = ["language", "command"]

    def __init__(self, elem):
        self.command = elem.text
        self.language = elem.get("language", "cheetah").lower()
        assert self.command
        assert self.language == "cheetah"


class MoveFileAction(BaseAction):
    action_type = "move_file"
    _keys = ["move_file"]

    def __init__(self, elem):
        self.source = elem.find("source").text
        self.destination = elem.find("destination").text

    def to_bash(self):
        return [f"mv {self.source} {self.destination}"], []


class MoveDirectoryFilesAction(BaseAction):
    action_type = "move_directory_files"
    _keys = ["source_directory", "destination_directory"]

    def __init__(self, elem):
        source_directory = elem.find("source_directory").text
        destination_directory = elem.find("destination_directory").text
        self.source_directory = source_directory
        self.destination_directory = destination_directory

    def to_bash(self):
        return [f"mv {self.source_directory}/* {self.destination_directory}/"], []


class SetEnvironmentAction(BaseAction):
    action_type = "set_environment"
    _keys = ["variables"]

    def __init__(self, elem):
        variables = []
        var_els = elem.findall("environment_variable")
        assert var_els is not None
        for ev_elem in var_els:
            var = SetVariable(ev_elem)
            variables.append(var)
        self.variables = variables
        assert self.variables

    def to_bash(self):
        answer = []
        for var in self.variables:
            # Expand $INSTALL_DIR here?
            if var.action == "set_to":
                answer.append(f"export {var.name}={var.raw_value}")
            elif var.action == "prepend_to":
                answer.append(f"export {var.name}={var.raw_value}:${var.name}")
            elif var.action == "append_to":
                answer.append(f"export {var.name}=${var.name}:{var.raw_value}")
            else:
                raise ValueError(f"Undefined environment variable action {var.action!r}")
        return answer, answer  # Actions needed in env.sh here!


class ChmodAction(BaseAction):
    action_type = "chmod"
    _keys = ["mods"]

    def __init__(self, elem):
        mods = []
        file_els = elem.findall("file")
        assert file_els is not None
        for mod_elem in file_els:
            mod = {}
            mod["mode"] = mod_elem.attrib["mode"]
            mod["target"] = mod_elem.text
            mods.append(mod)
        self.mods = mods
        assert self.mods

    def to_bash(self):
        return [f"chmod {m['mode']} {m['target']}" for m in self.mods], []


class MakeInstallAction(BaseAction):
    action_type = "make_install"
    _keys = []  # type: List[str]

    def __init__(self, elem):
        pass

    def to_bash(self):
        return ["make install"], []


class AutoconfAction(BaseAction):
    action_type = "autoconf"
    _keys = ["options"]

    def __init__(self, elem):
        self.options = elem.text

    def to_bash(self):
        if self.options:
            raise NotImplementedError("Options with action autoconf not implemented yet.")
        return ["./configure", "make", "make install"], []


class ChangeDirectoryAction(BaseAction):
    action_type = "change_directory"
    _keys = ["directory"]

    def __init__(self, elem):
        self.directory = elem.text
        assert self.directory

    def to_bash(self):
        return [f"cd {self.directory}"], []


class MakeDirectoryAction(BaseAction):
    action_type = "make_directory"
    _keys = ["directory"]

    def __init__(self, elem):
        self.directory = elem.text

    def to_bash(self):
        return [f"mkdir -p {self.directory}"], []


class SetupPerlEnvironmentAction(BaseAction):
    action_type = "setup_perl_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetupRubyEnvironmentAction(BaseAction):
    action_type = "setup_ruby_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetupPythonEnvironmentAction(BaseAction):
    action_type = "setup_python_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetupVirtualenvAction(BaseAction):
    action_type = "setup_virtualenv"
    _keys = ["use_requirements_file", "python", "requirements"]

    def __init__(self, elem):
        use_reqs = elem.attrib.get("use_requirements_file", "True")
        self.use_requirements_file = asbool(use_reqs)
        self.python = elem.get("python", "python")
        self.requirements = elem.text or "requirements.txt"


class SetupREnvironmentAction(BaseAction):
    action_type = "setup_r_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetEnvironmentForInstallAction(BaseAction):
    action_type = "set_environment_for_install"

    def __init__(self, elem):
        pass

    def to_bash(self):
        # TODO - How could we resolve/check the dependencies?
        return ['echo "WARNING: Assuming packages already installed!"'], []


class SetVariable:
    def __init__(self, elem):
        self.action = elem.attrib["action"]
        self.name = elem.attrib["name"]
        self.raw_value = elem.text


truthy = frozenset(["true", "yes", "on", "y", "t", "1"])
falsy = frozenset(["false", "no", "off", "n", "f", "0"])


def asbool(obj):
    if isinstance(obj, str):
        obj = obj.strip().lower()
        if obj in truthy:
            return True
        elif obj in falsy:
            return False
        else:
            raise ValueError(f"String is not true/false: {obj!r}")
    return bool(obj)


action_classes = [
    DownloadByUrlAction,
    DownloadFileAction,
    DownloadBinary,
    ShellCommandAction,
    TemplateShellCommandAction,
    MoveFileAction,
    MoveDirectoryFilesAction,
    SetEnvironmentAction,
    ChmodAction,
    MakeInstallAction,
    AutoconfAction,
    ChangeDirectoryAction,
    MakeDirectoryAction,
    SetupPerlEnvironmentAction,
    SetupRubyEnvironmentAction,
    SetupPythonEnvironmentAction,
    SetupVirtualenvAction,
    SetupREnvironmentAction,
    SetEnvironmentForInstallAction,
]

actions_by_type = dict(map(lambda c: (c.action_type, c), action_classes))
