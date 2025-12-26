import os
import unittest
from typing import (
    Any,
    Literal,
)

import bioblend
from bioblend.galaxy import GalaxyInstance
from . import test_util

bioblend.set_stream_logger("test", level="INFO")

BIOBLEND_TEST_JOB_TIMEOUT = int(os.environ.get("BIOBLEND_TEST_JOB_TIMEOUT", "60"))


@test_util.skip_unless_galaxy()
class GalaxyTestBase(unittest.TestCase):
    gi: GalaxyInstance

    @classmethod
    def setUpClass(cls) -> None:
        galaxy_key = os.environ["BIOBLEND_GALAXY_API_KEY"]
        galaxy_url = os.environ["BIOBLEND_GALAXY_URL"]
        cls.gi = GalaxyInstance(url=galaxy_url, key=galaxy_key)

    def _test_dataset(self, history_id: str, contents: str = "1\t2\t3", **kwargs: Any) -> str:
        tool_output = self.gi.tools.paste_content(contents, history_id, **kwargs)
        return tool_output["outputs"][0]["id"]

    def _wait_and_verify_dataset(
        self, dataset_id: str, expected_contents: bytes, timeout_seconds: float = BIOBLEND_TEST_JOB_TIMEOUT
    ) -> None:
        dataset_contents = self.gi.datasets.download_dataset(dataset_id, maxwait=timeout_seconds)
        assert dataset_contents == expected_contents

    def _run_random_lines1(
        self, history_id: str, dataset_id: str, input_format: Literal["21.01", "legacy"] = "legacy"
    ) -> dict[str, Any]:
        tool_inputs = {
            "num_lines": "1",
            "input": {"src": "hda", "id": dataset_id},
        }
        if input_format == "21.01":
            tool_inputs.update({"seed_source": {"seed_source_selector": "set_seed", "seed": "asdf"}})
        else:
            # legacy format
            tool_inputs.update({"seed_source|seed_source_selector": "set_seed", "seed_source|seed": "asdf"})
        return self.gi.tools.run_tool(
            history_id=history_id, tool_id="random_lines1", tool_inputs=tool_inputs, input_format=input_format
        )
