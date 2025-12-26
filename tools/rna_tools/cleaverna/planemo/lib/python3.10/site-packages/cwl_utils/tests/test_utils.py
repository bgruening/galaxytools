from urllib.parse import urlparse

from cwl_utils.utils import resolved_path


def test_resoled_path() -> None:
    base_url = urlparse(
        "schemas/bclconvert-run-configuration/2.0.0--4.0.3/bclconvert-run-configuration__2.0.0--4.0.3.yaml"
    )
    link = "../../../schemas/samplesheet/2.0.0--4.0.3/samplesheet__2.0.0--4.0.3.yaml#samplesheet"
    rpath = resolved_path(base_url, link)
    assert rpath == urlparse(
        "schemas/samplesheet/2.0.0--4.0.3/samplesheet__2.0.0--4.0.3.yaml#samplesheet"
    )
