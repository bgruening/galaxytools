# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def generate_instance_commands():
    """The Instance client holds the Singularity Instance command group
    The levels of verbosity (debug and quiet) are passed from the main
    client via the environment variable MESSAGELEVEL.

    """
    from spython.instance import Instance

    # run_command uses run_cmd, but wraps to catch error
    from spython.main.base.command import init_command, run_command
    from spython.main.base.generate import RobotNamer
    from spython.main.base.logger import println
    from spython.main.instances import list_instances

    from .logs import _logs, error_logs, output_logs
    from .start import start
    from .stop import stop

    Instance.RobotNamer = RobotNamer()
    Instance._init_command = init_command
    Instance.run_command = run_command
    Instance._list = list_instances  # list command is used to get metadata
    Instance._println = println
    Instance.start = start  # intended to be called on init, not by user
    Instance.stop = stop
    Instance.error_logs = error_logs
    Instance.output_logs = output_logs
    Instance._logs = _logs

    # Give an instance the ability to breed :)
    Instance.instance = Instance

    return Instance
