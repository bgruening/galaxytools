from . import GalaxyTestBase


class TestGalaxyToolData(GalaxyTestBase.GalaxyTestBase):
    def test_get_data_tables(self):
        tables = self.gi.tool_data.get_data_tables()
        for table in tables:
            assert table["name"] is not None

    def test_show_data_table(self):
        tables = self.gi.tool_data.get_data_tables()
        table = self.gi.tool_data.show_data_table(tables[0]["name"])
        assert table["columns"] is not None
        assert table["fields"] is not None
        assert table["name"] is not None
