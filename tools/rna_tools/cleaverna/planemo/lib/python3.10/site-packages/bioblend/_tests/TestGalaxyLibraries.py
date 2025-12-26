import os
import tempfile
from typing import Any

from . import (
    GalaxyTestBase,
    test_util,
)

FOO_DATA = "foo\nbar\n"


def listify(item: Any) -> list:
    # Slightly simplified version of listify() from
    # https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/util/__init__.py
    if not item:
        return []
    elif isinstance(item, (list, tuple)):
        return list(item)
    elif isinstance(item, str) and item.count(","):
        return [token.strip() for token in item.split(",")]
    return [item]


class TestGalaxyLibraries(GalaxyTestBase.GalaxyTestBase):
    def setUp(self):
        super().setUp()
        self.name = "automated test library"
        self.library = self.gi.libraries.create_library(
            self.name, description="automated test", synopsis="automated test synopsis"
        )

    def tearDown(self):
        self.gi.libraries.delete_library(self.library["id"])

    def test_create_library(self):
        assert self.library["name"] == self.name
        assert self.library["id"] is not None

    def test_get_libraries(self):
        libraries_with_name = self.gi.libraries.get_libraries(name=self.name)
        assert len([lib for lib in libraries_with_name if lib["id"] == self.library["id"]]) == 1

        deleted_name = "deleted test library"
        deleted_library_id = self.gi.libraries.create_library(
            deleted_name, description="a deleted library", synopsis="automated test synopsis"
        )["id"]
        self.gi.libraries.delete_library(deleted_library_id)
        deleted_libraries_with_name = self.gi.libraries.get_libraries(name=deleted_name, deleted=True)
        assert len([lib for lib in deleted_libraries_with_name if lib["id"] == deleted_library_id]) == 1

        all_non_deleted_libraries = self.gi.libraries.get_libraries(deleted=False)
        assert len([lib for lib in all_non_deleted_libraries if lib["id"] == self.library["id"]]) == 1
        assert [lib for lib in all_non_deleted_libraries if lib["id"] == deleted_library_id] == []

        all_deleted_libraries = self.gi.libraries.get_libraries(deleted=True)
        assert [lib for lib in all_deleted_libraries if lib["id"] == self.library["id"]] == []
        assert len([lib for lib in all_deleted_libraries if lib["id"] == deleted_library_id]) == 1

        all_libraries = self.gi.libraries.get_libraries(deleted=None)
        assert len([lib for lib in all_libraries if lib["id"] == self.library["id"]]) == 1
        assert len([lib for lib in all_libraries if lib["id"] == deleted_library_id]) == 1

    def test_show_library(self):
        library_data = self.gi.libraries.show_library(self.library["id"])
        assert self.library["id"] == library_data["id"]
        assert self.library["name"] == library_data["name"]

    def test_upload_file_from_url(self):
        url = "https://zenodo.org/record/582600/files/wildtype.fna?download=1"
        ret = self.gi.libraries.upload_file_from_url(self.library["id"], url)
        assert len(ret) == 1
        ldda_dict = ret[0]
        assert ldda_dict["name"] == url

    def test_upload_file_contents(self):
        ret = self.gi.libraries.upload_file_contents(self.library["id"], FOO_DATA)
        assert len(ret) == 1
        ldda_dict = ret[0]
        assert ldda_dict["name"] == "Pasted Entry"

    def test_upload_file_from_local_path(self):
        with tempfile.NamedTemporaryFile(mode="w", prefix="bioblend_test_") as f:
            f.write(FOO_DATA)
            f.flush()
            filename = f.name
            ret = self.gi.libraries.upload_file_from_local_path(self.library["id"], filename)
        assert len(ret) == 1
        ldda_dict = ret[0]
        assert ldda_dict["name"] == os.path.basename(filename)

    # def test_upload_file_from_server(self):
    #     pass

    def test_upload_from_galaxy_filesystem(self):
        bnames = [f"f{i}.txt" for i in range(2)]
        with tempfile.TemporaryDirectory() as tempdir:
            fnames = [os.path.join(tempdir, _) for _ in bnames]
            for fn in fnames:
                with open(fn, "w") as f:
                    f.write(FOO_DATA)
            filesystem_paths = "\n".join(fnames)
            ret = self.gi.libraries.upload_from_galaxy_filesystem(self.library["id"], filesystem_paths)
            for fn, dataset_dict in zip(fnames, ret):
                dataset = self.gi.libraries.wait_for_dataset(self.library["id"], dataset_dict["id"])
                assert dataset["state"] == "ok"
                assert dataset["name"] == os.path.basename(fn)
            ret = self.gi.libraries.upload_from_galaxy_filesystem(
                self.library["id"], filesystem_paths, link_data_only="link_to_files"
            )
            for fn, dataset_dict in zip(fnames, ret):
                dataset = self.gi.libraries.wait_for_dataset(self.library["id"], dataset_dict["id"])
                assert dataset["state"] == "ok"
                assert dataset["name"] == os.path.basename(fn)

    def test_copy_from_dataset(self):
        history = self.gi.histories.create_history()
        dataset_id = self._test_dataset(history["id"])
        self.gi.libraries.copy_from_dataset(self.library["id"], dataset_id, message="Copied from dataset")

    def test_update_dataset(self):
        library_id = self.library["id"]
        dataset1 = self.gi.libraries.upload_file_contents(library_id, FOO_DATA)
        updated_dataset = self.gi.libraries.update_library_dataset(
            dataset1[0]["id"], name="Modified name", misc_info="Modified the name succesfully"
        )
        assert updated_dataset["name"] == "Modified name"
        assert updated_dataset["misc_info"] == "Modified the name succesfully"

    def test_library_permissions(self):
        current_user = self.gi.users.get_current_user()
        user_id_list_new = [current_user["id"]]
        self.gi.libraries.set_library_permissions(
            self.library["id"],
            access_in=user_id_list_new,
            modify_in=user_id_list_new,
            add_in=user_id_list_new,
            manage_in=user_id_list_new,
        )
        ret = self.gi.libraries.get_library_permissions(self.library["id"])
        assert {_[1] for _ in ret["access_library_role_list"]} == set(user_id_list_new)
        assert {_[1] for _ in ret["modify_library_role_list"]} == set(user_id_list_new)
        assert {_[1] for _ in ret["add_library_item_role_list"]} == set(user_id_list_new)
        assert {_[1] for _ in ret["manage_library_role_list"]} == set(user_id_list_new)

    def test_dataset_permissions(self):
        current_user = self.gi.users.get_current_user()
        user_id_list_new = [current_user["id"]]
        library_id = self.library["id"]
        dataset1 = self.gi.libraries.upload_file_contents(library_id, FOO_DATA)
        ret = self.gi.libraries.set_dataset_permissions(
            dataset1[0]["id"], access_in=user_id_list_new, modify_in=user_id_list_new, manage_in=user_id_list_new
        )
        assert {_[1] for _ in ret["access_dataset_roles"]} == set(user_id_list_new)
        assert {_[1] for _ in ret["modify_item_roles"]} == set(user_id_list_new)
        assert {_[1] for _ in ret["manage_dataset_roles"]} == set(user_id_list_new)
        # test get_dataset_permissions
        ret_get = self.gi.libraries.get_dataset_permissions(dataset1[0]["id"])
        assert {_[1] for _ in ret_get["access_dataset_roles"]} == set(user_id_list_new)
        assert {_[1] for _ in ret_get["modify_item_roles"]} == set(user_id_list_new)
        assert {_[1] for _ in ret_get["manage_dataset_roles"]} == set(user_id_list_new)

    @test_util.skip_unless_galaxy("release_19.09")
    def test_upload_file_contents_with_tags(self):
        datasets = self.gi.libraries.upload_file_contents(self.library["id"], FOO_DATA, tags=["name:foobar", "barfoo"])
        dataset_show = self.gi.libraries.show_dataset(self.library["id"], datasets[0]["id"])
        assert listify(dataset_show["tags"]) == ["name:foobar", "barfoo"]

    @test_util.skip_unless_galaxy("release_19.09")
    def test_update_dataset_tags(self):
        datasets = self.gi.libraries.upload_file_contents(self.library["id"], FOO_DATA)
        dataset_show = self.gi.libraries.show_dataset(self.library["id"], datasets[0]["id"])
        assert listify(dataset_show["tags"]) == []

        updated_dataset = self.gi.libraries.update_library_dataset(datasets[0]["id"], tags=["name:foobar", "barfoo"])
        dataset_show = self.gi.libraries.show_dataset(self.library["id"], updated_dataset["id"])

        assert listify(dataset_show["tags"]) == ["name:foobar", "barfoo"]

    def test_folders(self):
        folder1_name = "folder 1"
        folder1 = self.gi.libraries.create_folder(self.library["id"], folder1_name)[0]
        assert folder1["name"] == folder1_name
        folder2_name = "folder 2"
        folder2 = self.gi.libraries.create_folder(self.library["id"], folder2_name)[0]
        assert folder2["name"] == folder2_name

        folders = self.gi.libraries.get_folders(self.library["id"])
        assert len(folders) == 3  # includes the root ("/") folder

        folder2_path = f"/{folder2_name}"
        folders = self.gi.libraries.get_folders(self.library["id"], name=folder2_path)
        assert len(folders) == 1
        assert folders[0]["name"] == folder2_path
