"""
A base representation of an instance of Galaxy
"""

from typing import Optional

from bioblend.galaxy import (
    config,
    container_resolution,
    dataset_collections,
    datasets,
    datatypes,
    folders,
    forms,
    ftpfiles,
    genomes,
    groups,
    histories,
    invocations,
    jobs,
    libraries,
    quotas,
    roles,
    tool_data,
    tool_dependencies,
    tools,
    toolshed,
    users,
    visual,
    workflows,
)
from bioblend.galaxyclient import GalaxyClient


class GalaxyInstance(GalaxyClient):
    def __init__(
        self,
        url: str,
        key: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None,
        *,
        token: Optional[str] = None,
        verify: bool = True,
        user_agent: Optional[str] = None,
    ) -> None:
        """
        A base representation of a connection to a Galaxy instance, identified
        by the server URL and user credentials.

        After you have created a ``GalaxyInstance`` object, access various
        modules via the class fields. For example, to work with histories and
        get a list of all the user's histories, the following should be done::

            from bioblend import galaxy

            gi = galaxy.GalaxyInstance(url='http://127.0.0.1:8000', key='your_api_key')

            hl = gi.histories.get_histories()

        :type url: str
        :param url: A FQDN or IP for a given instance of Galaxy. For example:
                    http://127.0.0.1:8080 . If a Galaxy instance is served under
                    a prefix (e.g., http://127.0.0.1:8080/galaxy/), supply the
                    entire URL including the prefix (note that the prefix must
                    end with a slash). If a Galaxy instance has HTTP Basic
                    authentication with username and password, then the
                    credentials should be included in the URL, e.g.
                    http://user:pass@host:port/galaxy/

        :type key: str
        :param key: User's API key for the given instance of Galaxy, obtained
                    from the user preferences. If a key is not supplied, an
                    email address and password must be and the key will
                    automatically be created for the user.

        :type email: str
        :param email: Galaxy e-mail address corresponding to the user.
                      Ignored if key is supplied directly.

        :type password: str
        :param password: Password of Galaxy account corresponding to the above
                         e-mail address. Ignored if key is supplied directly.

        :type token: str
        :param token: An OIDC access token obtained from an OIDC provider
                      configured in `oidc_backends_config.xml`. Can be used
                      as a substitue for an API key. You must make sure the access
                      token has not expired when making an API call.

        :param verify: Whether to verify the server's TLS certificate
        :type verify: bool
        """
        super().__init__(
            url, key=key, email=email, password=password, token=token, verify=verify, user_agent=user_agent
        )
        self.libraries = libraries.LibraryClient(self)
        self.histories = histories.HistoryClient(self)
        self.workflows = workflows.WorkflowClient(self)
        self.invocations = invocations.InvocationClient(self)
        self.datasets = datasets.DatasetClient(self)
        self.dataset_collections = dataset_collections.DatasetCollectionClient(self)
        self.users = users.UserClient(self)
        self.genomes = genomes.GenomeClient(self)
        self.tools = tools.ToolClient(self)
        self.toolshed = toolshed.ToolShedClient(self)
        self.toolShed = self.toolshed  # historical alias
        self.config = config.ConfigClient(self)
        self.container_resolution = container_resolution.ContainerResolutionClient(self)
        self.visual = visual.VisualClient(self)
        self.quotas = quotas.QuotaClient(self)
        self.groups = groups.GroupsClient(self)
        self.roles = roles.RolesClient(self)
        self.datatypes = datatypes.DatatypesClient(self)
        self.jobs = jobs.JobsClient(self)
        self.forms = forms.FormsClient(self)
        self.ftpfiles = ftpfiles.FTPFilesClient(self)
        self.tool_data = tool_data.ToolDataClient(self)
        self.folders = folders.FoldersClient(self)
        self.tool_dependencies = tool_dependencies.ToolDependenciesClient(self)

    def __repr__(self) -> str:
        """
        A nicer representation of this GalaxyInstance object
        """
        return f"GalaxyInstance object for Galaxy at {self.base_url}"
