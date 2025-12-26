"""
Test functions in bioblend.galaxy.tool_dependencies
"""

from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyToolDependencies(GalaxyTestBase.GalaxyTestBase):
    @test_util.skip_unless_galaxy("release_20.01")
    def test_summarize_toolbox(self):
        toolbox_summary = self.gi.tool_dependencies.summarize_toolbox()
        assert isinstance(toolbox_summary, list)
        assert len(toolbox_summary) > 0

        toolbox_summary_by_tool = self.gi.tool_dependencies.summarize_toolbox(index_by="tools")
        assert isinstance(toolbox_summary_by_tool, list)
        assert len(toolbox_summary_by_tool) > 0
        assert isinstance(toolbox_summary_by_tool[0], dict)
        assert "tool_ids" in toolbox_summary_by_tool[0]
        assert isinstance(toolbox_summary_by_tool[0]["tool_ids"], list)
        tool_id = toolbox_summary_by_tool[0]["tool_ids"][0]

        toolbox_summary_select_tool_ids = self.gi.tool_dependencies.summarize_toolbox(
            index_by="tools", tool_ids=[tool_id]
        )
        assert isinstance(toolbox_summary_select_tool_ids, list)
        assert len(toolbox_summary_select_tool_ids) == 1
        assert toolbox_summary_select_tool_ids[0]["tool_ids"][0] == tool_id

    @test_util.skip_unless_galaxy("release_20.01")
    def test_unused_dependency_paths(self):
        unused_paths = self.gi.tool_dependencies.unused_dependency_paths()
        assert isinstance(unused_paths, list)

    @test_util.skip_unless_galaxy("release_20.01")
    def test_delete_unused_dependency_paths(self):
        self.gi.tool_dependencies.delete_unused_dependency_paths(paths=[])
