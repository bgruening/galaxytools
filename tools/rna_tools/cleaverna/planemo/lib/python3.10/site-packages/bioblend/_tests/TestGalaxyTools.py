""" """

import os
from typing import Any

from bioblend.galaxy.tools.inputs import (
    conditional,
    dataset,
    inputs,
    repeat,
)
from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyTools(GalaxyTestBase.GalaxyTestBase):
    def test_get_tools(self):
        # Test requires target Galaxy is configured with at least one tool.
        tools = self.gi.tools.get_tools()
        assert len(tools) > 0
        assert all(map(self._assert_is_tool_rep, tools))

    def test_get_tool_panel(self):
        # Test requires target Galaxy is configured with at least one tool
        # section.
        tool_panel = self.gi.tools.get_tool_panel()
        sections = [s for s in tool_panel if "elems" in s]
        assert len(sections) > 0
        assert all(map(self._assert_is_tool_rep, sections[0]["elems"]))

    def _assert_is_tool_rep(self, data):
        assert data["model_class"].endswith("Tool")
        # Special tools like SetMetadataTool may have different model_class
        # than Tool - but they all seem to end in tool.

        for key in ["name", "id", "version"]:
            assert key in data
        return True

    def test_paste_content(self):
        history = self.gi.histories.create_history(name="test_paste_data history")
        paste_text = "line 1\nline 2\rline 3\r\nline 4"
        tool_output = self.gi.tools.paste_content(paste_text, history["id"])
        assert len(tool_output["outputs"]) == 1
        # All lines in the resulting dataset should end with "\n"
        expected_contents = ("\n".join(paste_text.splitlines()) + "\n").encode()
        self._wait_and_verify_dataset(tool_output["outputs"][0]["id"], expected_contents)
        # Same with space_to_tab=True
        tool_output = self.gi.tools.paste_content(paste_text, history["id"], space_to_tab=True)
        assert len(tool_output["outputs"]) == 1
        expected_contents = ("\n".join("\t".join(_.split()) for _ in paste_text.splitlines()) + "\n").encode()
        self._wait_and_verify_dataset(tool_output["outputs"][0]["id"], expected_contents)

    def test_upload_file(self):
        history = self.gi.histories.create_history(name="test_upload_file history")

        fn = test_util.get_abspath("test_util.py")
        file_name = "test1"
        tool_output = self.gi.tools.upload_file(
            fn,
            # First param could be a regular path also of course...
            history_id=history["id"],
            file_name=file_name,
            dbkey="?",
            file_type="txt",
        )
        self._wait_for_and_verify_upload(tool_output, file_name, fn, expected_dbkey="?")

    def test_upload_file_dbkey(self):
        history = self.gi.histories.create_history(name="test_upload_file history")
        fn = test_util.get_abspath("test_util.py")
        file_name = "test1"
        dbkey = "hg19"
        tool_output = self.gi.tools.upload_file(
            fn,
            history_id=history["id"],
            file_name=file_name,
            dbkey=dbkey,
            file_type="txt",
        )
        self._wait_for_and_verify_upload(tool_output, file_name, fn, expected_dbkey=dbkey)

    @test_util.skip_unless_tool("random_lines1")
    def test_run_random_lines(self):
        # Run second test case from randomlines.xml
        history_id = self.gi.histories.create_history(name="test_run_random_lines history")["id"]
        with open(test_util.get_abspath(os.path.join("data", "1.bed"))) as f:
            contents = f.read()
        dataset_id = self._test_dataset(history_id, contents=contents)
        tool_inputs = (
            inputs()
            .set("num_lines", "1")
            .set("input", dataset(dataset_id))
            .set("seed_source", conditional().set("seed_source_selector", "set_seed").set("seed", "asdf"))
        )
        tool_output = self.gi.tools.run_tool(history_id=history_id, tool_id="random_lines1", tool_inputs=tool_inputs)
        assert len(tool_output["outputs"]) == 1
        # TODO: Wait for results and verify has 1 line and is
        # chr5  131424298   131424460   CCDS4149.1_cds_0_0_chr5_131424299_f 0   +

    @test_util.skip_unless_tool("cat1")
    def test_run_cat1(self):
        history_id = self.gi.histories.create_history(name="test_run_cat1 history")["id"]
        dataset1_id = self._test_dataset(history_id, contents="1 2 3")
        dataset2_id = self._test_dataset(history_id, contents="4 5 6")
        dataset3_id = self._test_dataset(history_id, contents="7 8 9")
        tool_inputs = (
            inputs()
            .set("input1", dataset(dataset1_id))
            .set(
                "queries",
                repeat()
                .instance(inputs().set("input2", dataset(dataset2_id)))
                .instance(inputs().set("input2", dataset(dataset3_id))),
            )
        )
        tool_output = self.gi.tools.run_tool(history_id=history_id, tool_id="cat1", tool_inputs=tool_inputs)
        assert len(tool_output["outputs"]) == 1
        # TODO: Wait for results and verify it has 3 lines - 1 2 3, 4 5 6,
        # and 7 8 9.

    @test_util.skip_unless_tool("CONVERTER_fasta_to_bowtie_color_index")
    def test_tool_dependency_install(self):
        installed_dependencies = self.gi.tools.install_dependencies("CONVERTER_fasta_to_bowtie_color_index")
        assert any(
            True for d in installed_dependencies if d.get("name") == "bowtie" and d.get("dependency_type") == "conda"
        ), f"installed_dependencies is {installed_dependencies}"
        status = self.gi.tools.uninstall_dependencies("CONVERTER_fasta_to_bowtie_color_index")
        assert status[0]["model_class"] == "NullDependency", status

    @test_util.skip_unless_tool("CONVERTER_fasta_to_bowtie_color_index")
    def test_tool_requirements(self):
        tool_requirements = self.gi.tools.requirements("CONVERTER_fasta_to_bowtie_color_index")
        assert any(
            True
            for tr in tool_requirements
            if {"dependency_type", "version"} <= set(tr.keys()) and tr.get("name") == "bowtie"
        ), f"tool_requirements is {tool_requirements}"

    @test_util.skip_unless_tool("CONVERTER_fasta_to_bowtie_color_index")
    def test_reload(self):
        response = self.gi.tools.reload("CONVERTER_fasta_to_bowtie_color_index")
        assert isinstance(response, dict)
        assert "message" in response
        assert "id" in response["message"]

    @test_util.skip_unless_tool("sra_source")
    def test_get_citations(self):
        citations = self.gi.tools.get_citations("sra_source")
        assert len(citations) == 2

    def _wait_for_and_verify_upload(
        self, tool_output: dict[str, Any], file_name: str, fn: str, expected_dbkey: str = "?"
    ) -> None:
        assert len(tool_output["outputs"]) == 1
        output = tool_output["outputs"][0]
        assert output["name"] == file_name
        expected_contents = open(fn, "rb").read()
        self._wait_and_verify_dataset(output["id"], expected_contents)
        assert output["genome_build"] == expected_dbkey

    @test_util.skip_unless_tool("random_lines1")
    def test_get_tool_model(self):
        history_id = self.gi.histories.create_history(name="test_run_random_lines history")["id"]
        tool_model = self.gi.tools.build(tool_id="random_lines1", history_id=history_id)
        assert len(tool_model["inputs"]) == 3

    @test_util.skip_unless_tool("random_lines1")
    def test_get_tool_tests(self):
        tests = self.gi.tools.get_tool_tests("random_lines1")
        assert isinstance(tests, list)
        first_test = tests[0]
        assert first_test["tool_id"] == "random_lines1", first_test
        assert "inputs" in first_test
        assert "outputs" in first_test
