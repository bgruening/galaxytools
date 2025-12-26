"""

# Copyright (C) 2017-2022 Vanessa Sochat.

This Source Code Form is subject to the terms of the
Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""

import os
import pwd
import re
import shlex
import subprocess
import sys

from spython.logger import bot, decodeUtf8String


def _process_sudo_cmd(cmd, sudo, sudo_options):
    """
    Process the sudo command and honor adding environment (or not)
    """
    if sudo and sudo_options is not None:
        if isinstance(sudo_options, str):
            sudo_options = shlex.split(sudo_options)
        cmd = ["sudo"] + sudo_options + cmd
    elif sudo:
        cmd = ["sudo"] + cmd
    return [x for x in cmd if x]


def check_install(software="singularity", quiet=True):
    """
    check_install will attempt to run the singularity command, and
    return True if installed. The command line utils will not run
    without this check.
    """

    cmd = [software, "--version"]
    found = False

    try:
        version = run_command(cmd, quiet=True)
    except Exception:  # FileNotFoundError
        return found

    if version is not None:
        if version["return_code"] == 0:
            found = True

        if not quiet:
            version = version["message"]
            bot.info("Found %s version %s" % (software.upper(), version))

    return found


def get_singularity_version():
    """
    get the full singularity client version as reported by
    singularity --version [...]. For Singularity 3.x, this means:
    "singularity version 3.0.1-1"
    """
    version = os.environ.get("SPYTHON_SINGULARITY_VERSION", "")
    if version == "":
        try:
            version = run_command(["singularity", "--version"], quiet=True)
        except Exception:  # FileNotFoundError
            return version

        if version["return_code"] == 0:
            if version["message"]:
                version = version["message"][0].strip("\n")

    return version


def get_userhome():
    """
    Get the user home based on the effective uid
    """
    return pwd.getpwuid(os.getuid())[5]


def get_username():
    """
    Get the user name based on the effective uid
    """
    return pwd.getpwuid(os.getuid())[0]


def get_installdir():
    """
    Get_installdir returns the installation directory of the application
    """
    return os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


def stream_command(
    cmd,
    no_newline_regexp="Progess",
    sudo=False,
    sudo_options=None,
    output_type="stdout",
):
    """
    Stream a command (yield) back to the user, as each line is available.

    # Example usage:
    results = []
    for line in stream_command(cmd):
        print(line, end="")
        results.append(line)

    Parameters
    ==========
    cmd: the command to send, should be a list for subprocess
    no_newline_regexp: the regular expression to determine skipping a
                       newline. Defaults to finding Progress
    sudo_options: string or list of strings that will be passed as options to sudo

    """
    if output_type not in ["stdout", "stderr", "both"]:
        bot.exit(
            "Invalid output type %s. Must be stderr, stdout or both." % output_type
        )
    cmd = _process_sudo_cmd(cmd, sudo, sudo_options)

    stderr_pipe = subprocess.STDOUT if output_type == "both" else subprocess.PIPE

    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=stderr_pipe, universal_newlines=True
    )

    # Allow the runner to choose streaming output or error
    stream = process.stdout.readline
    if output_type == "stderr":
        stream = process.stderr.readline

    # Stream lines back to the caller
    for line in iter(stream, ""):
        if not re.search(no_newline_regexp, line):
            yield line

    # If there is an error, raise.
    process.stdout.close()
    return_code = process.wait()
    if return_code:
        # Some situations may return process without an attached stderr object
        # to read from
        if process.stderr:
            print(process.stderr.read(), file=sys.stderr)
        raise subprocess.CalledProcessError(return_code, cmd)


def run_command(
    cmd,
    sudo=False,
    capture=True,
    no_newline_regexp="Progess",
    quiet=False,
    sudo_options=None,
    environ=None,
    background=False,
):
    """
    run_command uses subprocess to send a command to the terminal. If
    capture is True, we use the parent stdout, so the progress bar (and
    other commands of interest) are piped to the user. This means we
    don't return the output to parse.

    Parameters
    ==========
    cmd: the command to send, should be a list for subprocess
    sudo: if needed, add to start of command
    no_newline_regexp: the regular expression to determine skipping a
                       newline. Defaults to finding Progress
    capture: if True, don't set stdout and have it go to console. This
             option can print a progress bar, but won't return the lines
             as output.
    sudo_options: string or list of strings that will be passed as options to sudo
    background: run in background and don't try to get output.
    """
    cmd = _process_sudo_cmd(cmd, sudo, sudo_options)

    stdout = None
    if capture:
        stdout = subprocess.PIPE

    # Use the parent stdout and stderr
    if background:
        subprocess.Popen(cmd, env=environ)
        return

    process = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=stdout, env=environ)

    lines = []
    found_match = False

    for line in process.communicate():
        if line:
            line = decodeUtf8String(line)
            lines.append(line)
            if re.search(no_newline_regexp, line) and found_match:
                if not quiet:
                    sys.stdout.write(line)
                found_match = True
            else:
                if not quiet:
                    sys.stdout.write(line)
                    print(line.rstrip())
                found_match = False

    output = {"message": lines, "return_code": process.returncode}

    return output


def format_container_name(name, special_characters=None):
    """
    format_container_name will take a name supplied by the user,
    remove all special characters (except for those defined by "special-characters"
    and return the new image name.
    """
    if special_characters is None:
        special_characters = []
    return "".join(e.lower() for e in name if e.isalnum() or e in special_characters)


def split_uri(container):
    """
    Split the uri of a container into the protocol and image part

    An empty protocol is returned if none found.
    A trailing slash is removed from the image part.
    """
    parts = container.split("://", 1)
    if len(parts) == 2:
        protocol, image = parts
    else:
        protocol = ""
        image = parts[0]
    return protocol, image.rstrip("/")


def remove_uri(container):
    """
    remove_uri will remove docker:// or shub:// or library:// from the uri
    """
    return split_uri(container)[1]
