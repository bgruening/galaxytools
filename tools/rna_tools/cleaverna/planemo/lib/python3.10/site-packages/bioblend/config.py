import configparser
import os
from typing import (
    IO,
    Optional,
)

BioBlendConfigPath = "/etc/bioblend.cfg"
BioBlendConfigLocations = [BioBlendConfigPath]
UserConfigPath = os.path.join(os.path.expanduser("~"), ".bioblend")
BioBlendConfigLocations.append(UserConfigPath)


class Config(configparser.ConfigParser):
    """
    BioBlend allows library-wide configuration to be set in external files.
    These configuration files can be used to specify access keys, for example.
    By default we use two locations for the BioBlend configurations:

    * System wide: ``/etc/bioblend.cfg``
    * Individual user: ``~/.bioblend`` (which works on both Windows and Unix)
    """

    def __init__(self, path: Optional[str] = None, fp: Optional[IO[str]] = None, do_load: bool = True) -> None:
        super().__init__({"working_dir": "/mnt/pyami", "debug": "0"})
        if do_load:
            if path:
                self.read([path])
            elif fp:
                self.read_file(fp)
            else:
                self.read(BioBlendConfigLocations)
