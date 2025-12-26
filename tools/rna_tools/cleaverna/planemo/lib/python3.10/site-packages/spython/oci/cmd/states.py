# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import json


def state(
    self, container_id=None, sudo=None, sync_socket=None, singularity_options=None
):
    """get the state of an OciImage, if it exists. The optional states that
    can be returned are created, running, stopped or (not existing).

    Equivalent command line example:
       singularity oci state <container_ID>

    Parameters
    ==========
    container_id: the id to get the state of.
    sudo: Add sudo to the command. If the container was created by root,
          you need sudo to interact and get its state.
    sync_socket: the path to the unix socket for state synchronization
    singularity_options: a list of options to provide to the singularity client

    Returns
    =======
    state: a parsed json of the container state, if exists. If the
           container is not found, None is returned.
    """
    sudo = self._get_sudo(sudo)
    container_id = self.get_container_id(container_id)

    # singularity oci state
    cmd = self._init_command("state", singularity_options)

    if sync_socket is not None:
        cmd = cmd + ["--sync-socket", sync_socket]

    # Finally, add the container_id
    cmd.append(container_id)

    # Get the instance state
    result = self._run_command(cmd, sudo=sudo, quiet=True)

    if result is not None:
        # If successful, a string is returned to parse
        if isinstance(result, str):
            return json.loads(result)


def _state_command(
    self, container_id=None, command="start", sudo=None, singularity_options=None
):
    """A generic state command to wrap pause, resume, kill, etc., where the
    only difference is the command. This function will be unwrapped if the
    child functions get more complicated (with additional arguments).

    Equivalent command line example:
       singularity oci <command> <container_ID>

    Parameters
    ==========
    container_id: the id to start.
    command: one of start, resume, pause, kill, defaults to start.
    singularity_options: a list of options to provide to the singularity client
    sudo: Add sudo to the command. If the container was created by root,
          you need sudo to interact and get its state.

    Returns
    =======
    return_code: the return code to indicate if the container was started.
    """
    sudo = self._get_sudo(sudo)
    container_id = self.get_container_id(container_id)

    # singularity oci state
    cmd = self._init_command(command, singularity_options)

    # Finally, add the container_id
    cmd.append(container_id)

    # Run the command, return return code
    return self._run_and_return(cmd, sudo)


def start(self, container_id=None, sudo=None, singularity_options=None):
    """start a previously invoked OciImage, if it exists.

    Equivalent command line example:
       singularity oci start <container_ID>

    Parameters
    ==========
    container_id: the id to start.
    sudo: Add sudo to the command. If the container was created by root,
          you need sudo to interact and get its state.

    Returns
    =======
    return_code: the return code to indicate if the container was started.
    """
    return self._state_command(
        container_id, sudo=sudo, singularity_options=singularity_options
    )


def kill(self, container_id=None, sudo=None, signal=None, singularity_options=None):
    """stop (kill) a started OciImage container, if it exists

    Equivalent command line example:
       singularity oci kill <container_ID>

    Parameters
    ==========
    container_id: the id to stop.
    signal: signal sent to the container (default SIGTERM)
    singularity_options: a list of options to provide to the singularity client
    sudo: Add sudo to the command. If the container was created by root,
          you need sudo to interact and get its state.

    Returns
    =======
    return_code: the return code to indicate if the container was killed.
    """
    sudo = self._get_sudo(sudo)
    container_id = self.get_container_id(container_id)

    # singularity oci state
    cmd = self._init_command("kill", singularity_options)

    # Finally, add the container_id
    cmd.append(container_id)

    # Add the signal, if defined
    if signal is not None:
        cmd = cmd + ["--signal", signal]

    # Run the command, return return code
    return self._run_and_return(cmd, sudo)


def resume(self, container_id=None, sudo=None, singularity_options=None):
    """resume a stopped OciImage container, if it exists

    Equivalent command line example:
       singularity oci resume <container_ID>

    Parameters
    ==========
    container_id: the id to stop.
    singularity_options: a list of options to provide to the singularity client
    sudo: Add sudo to the command. If the container was created by root,
          you need sudo to interact and get its state.

    Returns
    =======
    return_code: the return code to indicate if the container was resumed.
    """
    return self._state_command(
        container_id,
        command="resume",
        sudo=sudo,
        singularity_options=singularity_options,
    )


def pause(self, container_id=None, sudo=None, singularity_options=None):
    """pause a running OciImage container, if it exists

    Equivalent command line example:
       singularity oci pause <container_ID>

    Parameters
    ==========
    container_id: the id to stop.
    singularity_options: a list of options to provide to the singularity client
    sudo: Add sudo to the command. If the container was created by root,
          you need sudo to interact and get its state.

    Returns
    =======
    return_code: the return code to indicate if the container was paused.
    """
    return self._state_command(
        container_id,
        command="pause",
        sudo=sudo,
        singularity_options=singularity_options,
    )
