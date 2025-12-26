"""
Contains possible interactions with the Galaxy FTP Files
"""

from typing import TYPE_CHECKING

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class FTPFilesClient(Client):
    module = "ftp_files"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_ftp_files(self, deleted: bool = False) -> list[dict]:
        """
        Get a list of local files.

        :type deleted: bool
        :param deleted: Whether to include deleted files

        :rtype: list
        :return: A list of dicts with details on individual files on FTP
        """
        return self._get()
