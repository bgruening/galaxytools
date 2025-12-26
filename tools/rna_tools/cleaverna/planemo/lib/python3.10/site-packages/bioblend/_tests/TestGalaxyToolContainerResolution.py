"""
Test functions in bioblend.galaxy.container_resolution
"""

from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyContainerResolution(GalaxyTestBase.GalaxyTestBase):
    @test_util.skip_unless_galaxy("release_22.05")
    def test_get_container_resolvers(self):
        container_resolvers = self.gi.container_resolution.get_container_resolvers()
        assert isinstance(container_resolvers, list)
        assert len(container_resolvers) > 0
        assert isinstance(container_resolvers[0], dict)
        assert container_resolvers[0]["model_class"] == "ExplicitContainerResolver"
        assert container_resolvers[0]["resolver_type"] == "explicit"
        assert container_resolvers[0]["can_uninstall_dependencies"] is False
        assert container_resolvers[0]["builds_on_resolution"] is False

    @test_util.skip_unless_galaxy("release_22.05")
    def test_show_container_resolver(self):
        container_resolver = self.gi.container_resolution.show_container_resolver(0)
        print(container_resolver)
        assert isinstance(container_resolver, dict)
        assert container_resolver["model_class"] == "ExplicitContainerResolver"
        assert container_resolver["resolver_type"] == "explicit"
        assert container_resolver["can_uninstall_dependencies"] is False
        assert container_resolver["builds_on_resolution"] is False

    @test_util.skip_unless_galaxy("release_22.05")
    def test_resolve(self):
        tool = self.gi.container_resolution.resolve(tool_id="CONVERTER_parquet_to_csv")
        print(tool)
        assert isinstance(tool, dict)

        tool_requirements_only = self.gi.container_resolution.resolve(
            tool_id="CONVERTER_parquet_to_csv", requirements_only=True
        )
        assert isinstance(tool_requirements_only, dict)

    @test_util.skip_unless_galaxy("release_22.05")
    def test_resolve_toolbox(self):
        toolbox = self.gi.container_resolution.resolve_toolbox()
        assert isinstance(toolbox, list)
        assert len(toolbox) > 0
        assert isinstance(toolbox[0], dict)

        toolbox_by_tool_ids = self.gi.container_resolution.resolve_toolbox(tool_ids=[toolbox[0]["tool_id"]])
        assert isinstance(toolbox_by_tool_ids, list)
        assert len(toolbox_by_tool_ids) == 1
        assert isinstance(toolbox_by_tool_ids[0], dict)

        toolbox_by_resolver_type = self.gi.container_resolution.resolve_toolbox(resolver_type="mulled")
        assert isinstance(toolbox_by_resolver_type, list)
        assert len(toolbox_by_resolver_type) > 0
        assert isinstance(toolbox_by_resolver_type[0], dict)
        assert len(toolbox) == len(toolbox_by_resolver_type)
        for tool in toolbox_by_resolver_type:
            print(tool)
            assert (
                tool["status"]["dependency_type"] is None
                or tool["status"]["container_resolver"]["resolver_type"] == "mulled"
            )

        toolbox_by_container_type = self.gi.container_resolution.resolve_toolbox(container_type="docker")
        assert isinstance(toolbox_by_container_type, list)
        assert len(toolbox_by_container_type) > 0
        assert isinstance(toolbox_by_container_type[0], dict)
        assert len(toolbox) == len(toolbox_by_container_type)
        for tool in toolbox_by_container_type:
            assert tool["status"]["dependency_type"] is None or tool["status"]["dependency_type"] == "docker"
            assert (
                tool["status"]["dependency_type"] is None or tool["status"]["container_description"]["type"] == "docker"
            )

        toolbox_requirements_only = self.gi.container_resolution.resolve_toolbox(requirements_only=True)
        assert isinstance(toolbox_requirements_only, list)
        assert len(toolbox_requirements_only) > 0
        assert isinstance(toolbox_requirements_only[0], dict)
        assert len(toolbox) == len(toolbox_requirements_only)

        # TODO unless containers are available this may fallback to conda by default?
        #      depending on Galaxy's config
        # toolbox_by_index = self.gi.container_resolution.resolve_toolbox(tool_ids=[toolbox[0]['tool_id']], index=0, install=True)
        # assert isinstance(toolbox_by_index, list)
        # assert len(toolbox_by_index) > 0
        # assert isinstance(toolbox_by_index[0], dict)

    # TODO unless containers are available this may fallback to conda by default?
    #      depending on Galaxy's config
    # def test_resolve_toolbox_with_install(self):
    #     toolbox = self.gi.container_resolution.resolve_toolbox_with_install(tool_ids=[])
    #     assert isinstance(toolbox, list)
    #     assert len(toolbox) == 0
