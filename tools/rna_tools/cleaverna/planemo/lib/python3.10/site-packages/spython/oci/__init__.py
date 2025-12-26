# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


from spython.image import ImageBase
from spython.logger import bot


class OciImage(ImageBase):
    # Default functions of client don't use sudo
    sudo = False

    def __init__(
        self, container_id=None, bundle=None, create=True, sudo=True, **kwargs
    ):
        """An Oci Image is an Image Base with OCI functions appended

        Parameters
        ==========
        container_id: image uri to parse (required)
        bundle: a bundle directory to create a container from.
                the bundle should have a config.json at the root
        create: if the bundle is provided, create a container (default True)
        sudo: if init is called with or without sudo, keep a record and use
              for following commands unless sudo is provided to function.
        """
        super(OciImage, self).__init__()

        # Will typically be None, unless used outside of Client
        self.container_id = container_id
        self.protocol = "oci"
        self.sudo = sudo

        # If bundle is provided, create it
        if bundle is not None and container_id is not None and create:
            self.bundle = bundle
            self.create(bundle, container_id, **kwargs)

    # Unique resource identifier

    def get_container_id(self, container_id=None):
        """a helper function shared between functions that will return a
        container_id. First preference goes to a container_id provided by
        the user at runtime. Second preference goes to the container_id
        instantiated with the client.

        Parameters
        ==========
        container_id: image uri to parse (required)
        """

        # The user must provide a container_id, or have one with the client
        if container_id is None and self.container_id is None:
            bot.exit("You must provide a container_id.")

        # Choose whichever is not None, with preference for function provided
        container_id = container_id or self.container_id
        return container_id

    def get_uri(self):
        """return the image uri (oci://) along with it's name"""
        return self.__str__()

    # Naming

    def __str__(self):
        if self.container_id is not None:
            return "[singularity-python-oci:%s]" % self.container_id
        return "[singularity-python-oci]"

    def __repr__(self):
        return self.__str__()

    # Commands

    def _get_sudo(self, sudo=None):
        """if the client was initialized with sudo, remember this choice for
        later communication with the Oci Images. However, if the user provides
        a sudo argument (True or False) and not the default None, take
        preference to this argument.

        Parameters
        ==========
        sudo: if None, use self.sudo. Otherwise return sudo.
        """
        if sudo is None:
            sudo = self.sudo
        return sudo

    def _run_and_return(self, cmd, sudo=None):
        """Run a command, show the message to the user if quiet isn't set,
        and return the return code. This is a wrapper for the OCI client
        to run a command and easily return the return code value (what
        the user is ultimately interested in).

        Parameters
        ==========
        cmd: the command (list) to run.
        sudo: whether to add sudo or not.

        """
        sudo = self._get_sudo(sudo)
        result = self._run_command(cmd, sudo=sudo, quiet=True, return_result=True)

        # Successful return with no output
        if not result:
            return

        # Show the response to the user, only if not quiet.
        elif not self.quiet:
            bot.println(result["message"])

        # Return the state object to the user
        return result["return_code"]

    def _init_command(self, action, flags=None):
        """a wrapper to the base init_command, ensuring that "oci" is added
        to each command

        Parameters
        ==========
        action: the main action to perform (e.g., build)
        flags: one or more additional flags (e.g, volumes)
               not implemented yet.

        """
        from spython.main.base.command import init_command

        if not isinstance(action, list):
            action = [action]
        cmd = ["oci"] + action
        return init_command(self, cmd, flags)
