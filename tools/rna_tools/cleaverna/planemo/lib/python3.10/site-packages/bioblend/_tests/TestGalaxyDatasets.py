import shutil
import tempfile

import pytest

from bioblend import (
    ConnectionError,
    galaxy,
)
from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyDatasets(GalaxyTestBase.GalaxyTestBase):
    def setUp(self):
        super().setUp()
        self.history_id = self.gi.histories.create_history(name="TestDataset")["id"]
        self.dataset_contents = "line 1\nline 2\rline 3\r\nline 4"
        self.dataset_id = self._test_dataset(self.history_id, contents=self.dataset_contents)
        self.gi.datasets.wait_for_dataset(self.dataset_id)

    def tearDown(self):
        self.gi.histories.delete_history(self.history_id, purge=True)

    def test_show_nonexistent_dataset(self):
        with pytest.raises(ConnectionError):
            self.gi.datasets.show_dataset("nonexistent_id")

    def test_show_dataset(self):
        dataset = self.gi.datasets.show_dataset(self.dataset_id)
        assert dataset["id"] == self.dataset_id
        assert not dataset["deleted"]
        self.gi.histories.delete_dataset(self.history_id, self.dataset_id)
        dataset = self.gi.datasets.show_dataset(self.dataset_id)
        assert dataset["id"] == self.dataset_id
        assert dataset["deleted"]

    def test_download_dataset(self):
        with pytest.raises((TypeError, ConnectionError)):
            self.gi.datasets.download_dataset(None)  # type: ignore[call-overload]
        expected_contents = ("\n".join(self.dataset_contents.splitlines()) + "\n").encode()
        # download_dataset() with file_path=None is already tested in TestGalaxyTools.test_paste_content()
        # self._wait_and_verify_dataset(self.dataset_id, expected_contents)
        tempdir = tempfile.mkdtemp(prefix="bioblend_test_")
        try:
            downloaded_dataset = self.gi.datasets.download_dataset(
                self.dataset_id, file_path=tempdir, maxwait=GalaxyTestBase.BIOBLEND_TEST_JOB_TIMEOUT * 2
            )
            assert downloaded_dataset.startswith(tempdir)
            with open(downloaded_dataset, "rb") as f:
                assert f.read() == expected_contents
        finally:
            shutil.rmtree(tempdir)
        with tempfile.NamedTemporaryFile(prefix="bioblend_test_") as f:
            download_filename = self.gi.datasets.download_dataset(
                self.dataset_id,
                file_path=f.name,
                use_default_filename=False,
                maxwait=GalaxyTestBase.BIOBLEND_TEST_JOB_TIMEOUT,
            )
            assert download_filename == f.name
            f.flush()
            assert f.read() == expected_contents

    def test_get_datasets(self):
        datasets = self.gi.datasets.get_datasets()
        dataset_ids = [dataset["id"] for dataset in datasets]
        assert self.dataset_id in dataset_ids

    def test_get_datasets_history(self):
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id)
        assert len(datasets) == 1

    def test_get_datasets_limit_offset(self):
        datasets = self.gi.datasets.get_datasets(limit=1)
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, offset=1)
        assert datasets == []

    def test_get_datasets_name(self):
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, name="Pasted Entry")
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, name="Wrong Name")
        assert datasets == []

    @test_util.skip_unless_galaxy("release_20.05")
    def test_get_datasets_time(self):
        dataset = self.gi.datasets.show_dataset(self.dataset_id)
        ct = dataset["create_time"]
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, create_time_min=ct)
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, create_time_max=ct)
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, create_time_min="2100-01-01T00:00:00")
        assert datasets == []
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, create_time_max="2000-01-01T00:00:00")
        assert datasets == []

        ut = dataset["update_time"]
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, update_time_min=ut)
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, update_time_max=ut)
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, update_time_min="2100-01-01T00:00:00")
        assert datasets == []
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, update_time_max="2000-01-01T00:00:00")
        assert datasets == []

    @test_util.skip_unless_galaxy("release_20.05")
    def test_get_datasets_extension(self):
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, extension="txt")
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, extension="bam")
        assert datasets == []

    @test_util.skip_unless_galaxy("release_22.01")
    def test_get_datasets_extension_list(self):
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, extension=["bam", "txt"])
        assert len(datasets) == 1

    @test_util.skip_unless_galaxy("release_20.05")
    def test_get_datasets_state(self):
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, state="ok")
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, state="queued")
        assert datasets == []
        with pytest.raises(ConnectionError):
            self.gi.datasets.get_datasets(history_id=self.history_id, state="nonexistent_state")
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, state=["ok", "queued"])
        assert len(datasets) == 1

    @test_util.skip_unless_galaxy("release_20.05")
    def test_get_datasets_visible(self):
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, visible=True)
        assert len(datasets) == 1
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, visible=False)
        assert len(datasets) == 0

    def test_get_datasets_ordering(self):
        self.dataset_id2 = self._test_dataset(self.history_id, contents=self.dataset_contents)
        self.gi.datasets.wait_for_dataset(self.dataset_id2)
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, order="create_time-dsc")
        assert datasets[0]["id"] == self.dataset_id2
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, order="create_time-asc")
        assert datasets[0]["id"] == self.dataset_id
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, order="hid-dsc")
        assert datasets[0]["id"] == self.dataset_id2
        datasets = self.gi.datasets.get_datasets(history_id=self.history_id, order="hid-asc")
        assert datasets[0]["id"] == self.dataset_id

    def test_get_datasets_deleted(self):
        deleted_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, deleted=True)
        assert deleted_datasets == []
        self.gi.histories.delete_dataset(self.history_id, self.dataset_id)
        deleted_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, deleted=True)
        assert len(deleted_datasets) == 1
        purged_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, purged=True)
        assert purged_datasets == []
        self.gi.histories.delete_dataset(self.history_id, self.dataset_id, purge=True, wait=True)
        purged_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, purged=True)
        assert len(purged_datasets) == 1

    def test_get_datasets_tool_id_and_tag(self):
        cat1_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, tool_id="cat1")
        assert cat1_datasets == []
        upload1_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, tool_id="upload1")
        assert len(upload1_datasets) == 1
        self.gi.histories.update_dataset(self.history_id, self.dataset_id, tags=["test"])
        tagged_datasets = self.gi.datasets.get_datasets(history_id=self.history_id, tag="test")
        assert len(tagged_datasets) == 1

    def test_wait_for_dataset(self):
        history_id = self.gi.histories.create_history(name="TestWaitForDataset")["id"]
        dataset_contents = "line 1\nline 2\rline 3\r\nline 4"
        dataset_id = self._test_dataset(history_id, contents=dataset_contents)

        dataset = self.gi.datasets.wait_for_dataset(dataset_id)
        assert dataset["state"] == "ok"

        self.gi.histories.delete_history(history_id, purge=True)

    def test_dataset_permissions(self):
        anonymous_gi = galaxy.GalaxyInstance(url=self.gi.base_url, key=None)
        self.gi.datasets.publish_dataset(self.dataset_id, published=False)
        with pytest.raises(ConnectionError):
            anonymous_gi.datasets.show_dataset(self.dataset_id)
        self.gi.datasets.publish_dataset(self.dataset_id, published=True)
        # now dataset is public, i.e. accessible to anonymous users
        assert anonymous_gi.datasets.show_dataset(self.dataset_id)["id"] == self.dataset_id
        self.gi.datasets.publish_dataset(self.dataset_id, published=False)

        admin_user_id = self.gi.users.get_current_user()["id"]
        new_user, user_gi = test_util.new_user_gi(self.gi)
        sharing_role = self.gi.roles.create_role("sharing_role", "sharing_role", [new_user["id"], admin_user_id])["id"]
        with pytest.raises(ConnectionError):
            user_gi.datasets.show_dataset(self.dataset_id)
        self.gi.datasets.update_permissions(self.dataset_id, access_ids=[sharing_role], manage_ids=[sharing_role])
        assert user_gi.datasets.show_dataset(self.dataset_id)["id"] == self.dataset_id
        # anonymous access now fails because sharing is only with the shared user role
        with pytest.raises(ConnectionError):
            anonymous_gi.datasets.show_dataset(self.dataset_id)
