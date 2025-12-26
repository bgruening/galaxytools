# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os

from spython.logger import bot


def export(
    self,
    image_path,
    pipe=False,
    output_file=None,
    command=None,
    sudo=False,
    singularity_options=None,
):
    """export will export an image, sudo must be used. Since we have Singularity
    versions after 3, export is replaced with building into a sandbox.

    Parameters
    ==========
    image_path: full path to image
    pipe: export to pipe and not file (default, False)
    singularity_options: a list of options to provide to the singularity client
    output_file: if pipe=False, export tar to this file. If not specified,
    will generate temporary directory.
    """
    from spython.utils import check_install

    check_install()

    # If export is deprecated, we run a build
    bot.warning(
        "Export is not supported for Singularity 3.x. Building to sandbox instead."
    )

    if output_file is None:
        basename, _ = os.path.splitext(image_path)
        output_file = self._get_filename(basename, "sandbox", pwd=False)

    return self.build(
        recipe=image_path,
        image=output_file,
        sandbox=True,
        force=True,
        sudo=sudo,
        singularity_options=singularity_options,
    )
