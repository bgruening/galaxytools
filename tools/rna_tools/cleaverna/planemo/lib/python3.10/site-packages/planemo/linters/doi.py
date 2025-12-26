"""Tool linting module that lints Galaxy tools for their DOIs (if a DOI type citation is present)"""

import planemo.lint


def lint_tool_dois(tool_xml, lint_ctx):
    planemo.lint.lint_dois(tool_xml, lint_ctx)
