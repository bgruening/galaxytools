import os
import unittest

import bioblend
import bioblend.toolshed
from . import test_util


@test_util.skip_unless_toolshed()
class TestToolshed(unittest.TestCase):
    def setUp(self):
        toolshed_url = os.environ["BIOBLEND_TOOLSHED_URL"]
        self.ts = bioblend.toolshed.ToolShedInstance(url=toolshed_url)

    def test_categories_client(self):
        # get_categories
        categories = self.ts.categories.get_categories()
        assert "Assembly" in [c["name"] for c in categories]
        # we cannot test get_categories with deleted=True as it requires administrator status

        # show_category
        visualization_category_id = [c for c in categories if c["name"] == "Visualization"][0]["id"]
        visualization_category = self.ts.categories.show_category(visualization_category_id)
        assert visualization_category["description"] == "Tools for visualizing data"

        # get_repositories
        repositories = self.ts.categories.get_repositories(visualization_category_id)
        repositories_reversed = self.ts.categories.get_repositories(visualization_category_id, sort_order="desc")
        assert len(repositories["repositories"]) > 200
        assert {
            "deprecated",
            "description",
            "homepage_url",
            "id",
            "name",
            "owner",
            "remote_repository_url",
            "type",
            "update_time",
        } <= set(repositories["repositories"][0].keys())
        assert repositories["repositories"][0] == repositories_reversed["repositories"][-1]

    def test_repositories_client(self):
        # get_repositories
        repositories = self.ts.repositories.get_repositories()
        assert len(repositories) > 5000
        repository0 = repositories[0]
        for key in ("id", "name", "owner", "type", "description", "deprecated"):
            assert key in repository0

        repositories = self.ts.repositories.get_repositories(name="bam_to_sam", owner="devteam")
        assert len(repositories) == 1
        bam_to_sam_repo = repositories[0]
        assert bam_to_sam_repo["name"] == "bam_to_sam"
        assert bam_to_sam_repo["owner"] == "devteam"
        assert bam_to_sam_repo["type"] == "unrestricted"
        assert not bam_to_sam_repo["deprecated"]

        # search_repositories
        samtools_search = self.ts.repositories.search_repositories("samtools", page_size=5)
        assert int(samtools_search["total_results"]) > 20
        assert len(samtools_search["hits"]) == 5

        # show_repository
        show_bam_to_sam_repo = self.ts.repositories.show_repository(bam_to_sam_repo["id"])
        assert "SAM" in show_bam_to_sam_repo["long_description"]

        # test_create_repository
        # need to provide an API key to test this

        # test_update_repository
        # need to provide an API key to test this

    def test_repositories_revisions(self):
        # get_ordered_installable_revisions
        bam_to_sam_revisions = self.ts.repositories.get_ordered_installable_revisions("bam_to_sam", "devteam")
        assert len(bam_to_sam_revisions) >= 4

        # get_repository_revision_install_info
        bam_to_sam_revision_install_info = self.ts.repositories.get_repository_revision_install_info(
            "bam_to_sam", "devteam", bam_to_sam_revisions[0]
        )
        assert len(bam_to_sam_revision_install_info) == 3
        assert bam_to_sam_revision_install_info[0].get("model_class") == "Repository"
        assert bam_to_sam_revision_install_info[1].get("model_class") == "RepositoryMetadata"
        assert bam_to_sam_revision_install_info[2].get("model_class") is None

    def test_tools_client(self):
        # search_tools
        samtools_search = self.ts.tools.search_tools("samtools", page_size=5)
        assert int(samtools_search["page"]) == 1
        assert len(samtools_search["hits"]) == 5
        hit0_tool = samtools_search["hits"][0]["tool"]
        for key in ("id", "repo_owner_username", "repo_name", "name", "description"):
            assert key in hit0_tool
