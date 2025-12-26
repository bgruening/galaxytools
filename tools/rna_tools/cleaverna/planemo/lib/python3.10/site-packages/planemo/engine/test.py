from planemo.engine import engine_context
from planemo.galaxy.config import _find_test_data
from planemo.galaxy.test import handle_reports_and_summary
from planemo.runnable import for_paths


def test_runnables(ctx, runnables, original_paths=None, **kwds):
    """Return exit code indicating test or failure."""
    if kwds.get("update_test_data"):
        non_copied_runnables = for_paths(original_paths)
        kwds["test_data_target_dir"] = _find_test_data(non_copied_runnables, **kwds)
    with engine_context(ctx, **kwds) as engine:
        test_data = engine.test(runnables, test_timeout=kwds.get("test_timeout"))
        ctx.vlog(f"engine.test returning [{test_data}]")
        return handle_reports_and_summary(ctx, test_data.structured_data, kwds=kwds)
