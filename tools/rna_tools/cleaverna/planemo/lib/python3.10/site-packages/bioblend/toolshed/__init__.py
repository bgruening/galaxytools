"""
A base representation of an instance of Tool Shed
"""

from typing import Optional

from bioblend.galaxyclient import GalaxyClient
from bioblend.toolshed import (
    categories,
    repositories,
    tools,
)


class ToolShedInstance(GalaxyClient):
    def __init__(
        self,
        url: str,
        key: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None,
        *,
        verify: bool = True,
        user_agent: Optional[str] = None,
    ) -> None:
        """
        A base representation of a connection to a ToolShed instance, identified
        by the ToolShed URL and user credentials.

        After you have created a ``ToolShedInstance`` object, access various
        modules via the class fields. For example, to work with repositories and
        get a list of all public repositories, the following should be done::

            from bioblend import toolshed

            ts = toolshed.ToolShedInstance(url='https://testtoolshed.g2.bx.psu.edu')

            rl = ts.repositories.get_repositories()

            tools = ts.tools.search_tools('fastq')

        :type url: str
        :param url: A FQDN or IP for a given instance of ToolShed. For example:
                    https://testtoolshed.g2.bx.psu.edu . If a ToolShed instance
                    is served under a prefix (e.g.
                    http://127.0.0.1:8080/toolshed/), supply the entire URL
                    including the prefix (note that the prefix must end with a
                    slash).

        :type key: str
        :param key: If required, user's API key for the given instance of ToolShed,
                    obtained from the user preferences.

        :type email: str
        :param email: ToolShed e-mail address corresponding to the user.
                      Ignored if key is supplied directly.

        :type password: str
        :param password: Password of ToolShed account corresponding to the above
                         e-mail address. Ignored if key is supplied directly.

        :param verify: Whether to verify the server's TLS certificate
        :type verify: bool
        """
        super().__init__(url, key=key, email=email, password=password, verify=verify, user_agent=user_agent)
        self.categories = categories.ToolShedCategoryClient(self)
        self.repositories = repositories.ToolShedRepositoryClient(self)
        self.tools = tools.ToolShedToolClient(self)
