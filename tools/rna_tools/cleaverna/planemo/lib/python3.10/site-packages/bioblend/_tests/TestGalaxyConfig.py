from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyConfig(GalaxyTestBase.GalaxyTestBase):
    def test_get_config(self):
        response = self.gi.config.get_config()
        assert isinstance(response, dict)
        assert "brand" in response.keys()

    def test_get_version(self):
        response = self.gi.config.get_version()
        assert isinstance(response, dict)
        assert "version_major" in response.keys()

    def test_whoami(self):
        response = self.gi.config.whoami()
        assert isinstance(response, dict)
        assert "username" in response.keys()

    def test_reload_toolbox(self):
        response = self.gi.config.reload_toolbox()
        assert response is None

    @test_util.skip_unless_galaxy("release_24.0")
    def test_encode_decode_id(self):
        int_id = 42
        encoded_id = self.gi.config.encode_id(int_id)
        assert isinstance(encoded_id, str)
        decoded_id = self.gi.config.decode_id(encoded_id)
        assert isinstance(decoded_id, int)
        assert decoded_id == int_id
