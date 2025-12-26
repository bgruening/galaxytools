import uuid

from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyQuotas(GalaxyTestBase.GalaxyTestBase):
    def setUp(self):
        super().setUp()
        # Quota names must be unique, and they're impossible to delete
        # without accessing the database.
        self.quota_name = f"BioBlend-Test-Quota-{uuid.uuid4().hex}"
        self.quota = self.gi.quotas.create_quota(self.quota_name, "testing", "100 GB", "=", default="registered")

    def tearDown(self):
        self.gi.quotas.update_quota(self.quota["id"], default="registered")
        self.gi.quotas.update_quota(self.quota["id"], default="no")
        self.gi.quotas.delete_quota(self.quota["id"])

    def test_create_quota(self):
        quota = self.gi.quotas.show_quota(self.quota["id"])
        assert quota["id"] == self.quota["id"]
        assert quota["name"] == self.quota_name
        assert quota["bytes"] == 107374182400
        assert quota["operation"] == "="
        assert quota["description"] == "testing"

    def test_get_quotas(self):
        quotas = self.gi.quotas.get_quotas()
        assert self.quota["id"] in [quota["id"] for quota in quotas]

    def test_update_quota(self):
        response = self.gi.quotas.update_quota(
            self.quota["id"],
            name=self.quota_name + "-new",
            description="asdf",
            default="registered",
            operation="-",
            amount=".01 TB",
        )
        assert f"""Quota '{self.quota_name}' has been renamed to '{self.quota_name}-new'""" in response

        quota = self.gi.quotas.show_quota(self.quota["id"])
        assert quota["id"] == self.quota["id"]
        assert quota["name"] == self.quota_name + "-new"
        assert quota["bytes"] == 10995116277
        assert quota["operation"] == "-"
        assert quota["description"] == "asdf"

    def test_delete_undelete_quota(self):
        self.gi.quotas.update_quota(self.quota["id"], default="no")
        response = self.gi.quotas.delete_quota(self.quota["id"])
        assert response == "Deleted 1 quotas: " + self.quota_name
        response = self.gi.quotas.undelete_quota(self.quota["id"])
        assert response == "Undeleted 1 quotas: " + self.quota_name

    @test_util.skip_unless_galaxy("release_19.09")  # for user purging
    def test_update_non_default_quota(self):
        """
        Test updating a non default quota.
        Needs to use `default=None` (which is the default), `default="no"` will fail.
        """
        if self.gi.config.get_config()["use_remote_user"]:
            self.skipTest("This Galaxy instance is not configured to use local users")
        new_username = test_util.random_string()
        new_user_email = f"{new_username}@example.org"
        password = test_util.random_string(20)
        new_user = self.gi.users.create_local_user(new_username, new_user_email, password)

        quota = self.gi.quotas.create_quota(
            name="non_default_quota",
            description="testing",
            amount="100 GB",
            operation="+",
            in_users=[new_user["id"]],
        )
        self.gi.quotas.update_quota(quota["id"], amount="200 GB")

        if self.gi.config.get_config()["allow_user_deletion"]:
            self.gi.users.delete_user(new_user["id"])
            self.gi.users.delete_user(new_user["id"], purge=True)
