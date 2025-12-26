"""
Contains possible interactions with the Galaxy Datatype
"""

from typing import TYPE_CHECKING

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class DatatypesClient(Client):
    module = "datatypes"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_datatypes(self, extension_only: bool = False, upload_only: bool = False) -> list[str]:
        """
        Get the list of all installed datatypes.

        :type extension_only: bool
        :param extension_only: Return only the extension rather than the datatype name

        :type upload_only: bool
        :param upload_only: Whether to return only datatypes which can be uploaded

        :rtype: list
        :return: A list of datatype names.
          For example::

            ['snpmatrix',
             'snptest',
             'tabular',
             'taxonomy',
             'twobit',
             'txt',
             'vcf',
             'wig',
             'xgmml',
             'xml']
        """

        params: dict[str, bool] = {}
        if extension_only:
            params["extension_only"] = True

        if upload_only:
            params["upload_only"] = True

        return self._get(params=params)

    def get_sniffers(self) -> list[str]:
        """
        Get the list of all installed sniffers.

        :rtype: list
        :return: A list of sniffer names.
          For example::

            ['galaxy.datatypes.tabular:Vcf',
             'galaxy.datatypes.binary:TwoBit',
             'galaxy.datatypes.binary:Bam',
             'galaxy.datatypes.binary:Sff',
             'galaxy.datatypes.xml:Phyloxml',
             'galaxy.datatypes.xml:GenericXml',
             'galaxy.datatypes.sequence:Maf',
             'galaxy.datatypes.sequence:Lav',
             'galaxy.datatypes.sequence:csFasta']
        """
        url = self._make_url() + "/sniffers"
        return self._get(url=url)
