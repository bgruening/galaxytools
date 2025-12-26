# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def main(args, options, parser):
    # If we have options, first is image
    image = None
    if options:
        image = options.pop(0)

    lookup = {"ipython": ipython, "python": python, "bpython": run_bpython}

    shells = ["ipython", "python", "bpython"]

    # Otherwise present order of liklihood to have on system
    for shell in shells:
        try:
            return lookup[shell](image)
        except ImportError:
            pass


def prepare_client(image):
    """prepare a client to embed in a shell with recipe parsers and writers."""
    # The client will announce itself (backend/database) unless it's get
    from spython.main import get_client
    from spython.main.parse import parsers, writers

    client = get_client()

    if image:
        client.load(image)

    # Add recipe parsers
    client.parsers = parsers
    client.writers = writers
    return client


def ipython(image):
    """give the user an ipython shell"""
    client = prepare_client(image)  # noqa

    try:
        from IPython import embed
    except ImportError:
        return python(image)

    embed()


def run_bpython(image):
    """give the user a bpython shell"""
    client = prepare_client(image)

    try:
        import bpython
    except ImportError:
        return python(image)

    bpython.embed(locals_={"client": client})


def python(image):
    """give the user a python shell"""
    import code

    client = prepare_client(image)
    code.interact(local={"client": client})
