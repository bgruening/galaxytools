# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


from spython.logger import bot


def stop(
    self,
    name=None,
    sudo=False,
    sudo_options=None,
    timeout=None,
    singularity_options=None,
    quiet=True,
):
    """stop an instance. This is done by default when an instance is created.

    Parameters
    ==========
    name: a name for the instance
    sudo: if the user wants to run the command with sudo
    singularity_options: a list of options to provide to the singularity client
    quiet: Do not print verbose output.
    timeout: forcebly kill non-stopped instance after the
             timeout specified in seconds

    USAGE:
    singularity [...] instance.stop [...] <instance name>

    """
    from spython.utils import check_install, run_command

    check_install()

    subgroup = ["instance", "stop"]
    if timeout:
        subgroup += ["-t", str(timeout)]

    cmd = self._init_command(subgroup, singularity_options)

    # If name is provided assume referencing an instance
    instance_name = self.name
    if name is not None:
        instance_name = name
    cmd = cmd + [instance_name]

    # Print verbose output
    if not (quiet or self.quiet):
        bot.info(" ".join(cmd))

    output = run_command(cmd, sudo=sudo, sudo_options=sudo_options, quiet=True)

    if output["return_code"] != 0:
        message = "%s : return code %s" % (output["message"], output["return_code"])
        bot.error(message)
        return output["return_code"]

    return output["return_code"]
