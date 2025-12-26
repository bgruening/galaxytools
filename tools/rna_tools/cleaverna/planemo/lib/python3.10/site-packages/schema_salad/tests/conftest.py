import pytest

import schema_salad.schema


@pytest.fixture(autouse=True, scope="function")
def isolated_cache() -> None:
    """
    Clear the schema_salad metaschema cache.

    Auto-loaded (see autouse) fixture, loaded per test (function scope).
    Prevents issues when running multiple tests that load metaschemas
    multiple times or in parallel (`pytest-parallel`, `pytest-xdist`, etc).
    """
    schema_salad.schema.cached_metaschema = None
