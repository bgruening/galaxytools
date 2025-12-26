"""Tool linting module that lints Galaxy tools for their URLs"""

import planemo.lint


def lint_tool_urls(tool_source, lint_ctx):
    planemo.lint.lint_urls(tool_source.root, lint_ctx)
