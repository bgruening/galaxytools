# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def apps(self, image=None, full_path=False, root=""):
    """
    return list of SCIF apps in image. The Singularity software serves
    a scientific filesystem integration that will install apps to
    /scif/apps and associated data to /scif/data. For more information
    about SCIF, see https://sci-f.github.io. Note that this seems
    to be deprecated in Singularity 3.x.

    Parameters
    ==========
    full_path: if True, return relative to scif base folder
    image_path: full path to the image

    """
    from spython.utils import check_install

    check_install()

    # No image provided, default to use the client's loaded image
    if image is None:
        image = self._get_uri()

    cmd = self._init_command("apps") + [image]
    output = self._run_command(cmd)

    if full_path:
        root = "/scif/apps/"

    if output:
        output = "".join(output).split("\n")
        output = ["%s%s" % (root, x) for x in output if x]

    return output
