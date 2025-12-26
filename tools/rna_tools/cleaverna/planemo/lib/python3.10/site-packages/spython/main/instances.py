# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json

from spython.logger import bot
from spython.utils import run_command


def list_instances(
    self,
    name=None,
    return_json=False,
    quiet=False,
    sudo=False,
    sudo_options=None,
    singularity_options=None,
):
    """
    List instances. For Singularity, this is provided as a command sub
    group.

    singularity instance list

    Return codes provided are different from standard linux:
    see https://github.com/singularityware/singularity/issues/1706
    Since we expect json output, we don't support older versions of Singularity.

    Parameters
    ==========
    return_json: return a json list of instances instead of objects (False)
    name: if defined, return the list for just one instance (used to ged pid)
    singularity_options: a list of options to provide to the singularity client

    Return Code  --   Reason
    0 -- Instances Found
    1 -- No Instances, libexecdir value not found, functions file not found
    255 -- Couldn't get UID

    """
    from spython.utils import check_install

    check_install()

    subgroup = ["instance", "list", "--json"]
    cmd = self._init_command(subgroup, singularity_options)

    # If the user has provided a name, we want to see a particular instance
    if name is not None:
        cmd.append(name)

    # Does the user want to see the command printed?
    if not (quiet or self.quiet):
        bot.info(" ".join(cmd))

    output = run_command(cmd, quiet=True, sudo=sudo, sudo_options=sudo_options)
    instances = []

    # Success, we have instances

    if output["return_code"] == 0:
        instances = json.loads(output["message"][0]).get("instances", {})

        # Does the user want instance objects instead?
        listing = []

        if not return_json:
            for i in instances:
                # If the user has provided a name, only add instance matches
                if name is not None:
                    if name != i["instance"]:
                        continue

                # Otherwise, add instances to the listing
                new_instance = self.instance(
                    pid=i.get("pid"),
                    ip_address=i.get("ip"),
                    name=i.get("instance") or i.get("daemon_name"),
                    log_err_path=i.get("logErrPath"),
                    log_out_path=i.get("logOutPath"),
                    image=i.get("img") or i.get("container_image"),
                    start=False,
                )

                listing.append(new_instance)
            instances = listing

    # Couldn't get UID

    elif output["return_code"] == 255:
        bot.error("Couldn't get UID")

    # Return code of 0
    else:
        bot.info("No instances found.")

    # If we are given a name, return just one
    if name is not None and instances and len(instances) == 1:
        instances = instances[0]

    return instances


def stopall(self, sudo=False, quiet=True, singularity_options=None):
    """
    Stop ALL instances. This command is only added to the command group
    as it doesn't make sense to call from a single instance

    Parameters
    ==========
    sudo: if the command should be done with sudo (exposes different set of
          instances)

    """
    from spython.utils import check_install

    check_install()

    subgroup = "instance.stop"

    if "version 3" in self.version():
        subgroup = ["instance", "stop"]

    cmd = self._init_command(subgroup, singularity_options)
    cmd = cmd + ["--all"]

    # Does the user want to see the command printed?
    if not (quiet or self.quiet):
        bot.info(" ".join(cmd))

    output = run_command(cmd, sudo=sudo, quiet=quiet)

    if output["return_code"] != 0:
        message = "%s : return code %s" % (output["message"], output["return_code"])
        bot.error(message)
        return output["return_code"]

    return output["return_code"]
