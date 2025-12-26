"""Utilities to help linting various targets."""

import os
import re
from typing import (
    Any,
    Dict,
    TYPE_CHECKING,
)
from urllib.request import urlopen

import requests
from galaxy.tool_util.lint import (
    LintContext,
    Linter,
)

from planemo.io import error
from planemo.shed import find_urls_for_xml
from planemo.xml import validation

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

REQUEST_TIMEOUT = 5


def build_lint_args(ctx: "PlanemoCliContext", **kwds) -> Dict[str, Any]:
    """Handle common report, error, and skip linting arguments."""
    report_level = kwds.get("report_level", "all")
    fail_level = kwds.get("fail_level", "warn")
    skip = kwds.get("skip", ctx.global_config.get("lint_skip"))
    if skip is None:
        skip = []
    if isinstance(skip, list):
        skip_types = skip
    else:
        skip_types = [s.strip() for s in skip.split(",")]

    for skip_file in kwds.get("skip_file", []):
        with open(skip_file) as f:
            for line in f.readlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                skip_types.append(line)

    linters = Linter.list_linters()
    linters.extend(["version_bumped"])
    invalid_skip_types = list(set(skip_types) - set(linters))
    if len(invalid_skip_types):
        error(f"Unknown linter type(s) {invalid_skip_types} in list of linters to be skipped. Known linters {linters}")

    lint_args: Dict[str, Any] = dict(
        level=report_level,
        fail_level=fail_level,
        skip_types=skip_types,
    )
    return lint_args


def setup_lint(ctx, **kwds):
    """Prepare lint_args and lint_ctx to begin linting a target."""
    lint_args = kwds.get("lint_args", None) or build_lint_args(ctx, **kwds)
    lint_ctx = LintContext(level=lint_args["level"], skip_types=lint_args["skip_types"])
    return lint_args, lint_ctx


def handle_lint_complete(lint_ctx, lint_args, failed=False):
    """Complete linting of a target and decide exit code."""
    if not failed:
        failed = lint_ctx.failed(lint_args["fail_level"])
    if failed:
        error("Failed linting")
    return 1 if failed else 0


def lint_dois(tool_xml, lint_ctx):
    """Find referenced DOIs and check they have valid with https://doi.org."""
    dois = find_dois_for_xml(tool_xml)
    for publication in dois:
        is_doi(publication, lint_ctx)


def find_dois_for_xml(tool_xml):
    dois = []
    for element in tool_xml.getroot().findall("citations"):
        for citation in list(element):
            if citation.tag == "citation" and citation.attrib.get("type", "") == "doi":
                dois.append(citation.text)
    return dois


def is_doi(publication_id, lint_ctx):
    """Check if dx.doi knows about the ``publication_id``."""
    base_url = "https://doi.org"
    if publication_id is None:
        lint_ctx.error("Empty DOI citation")
        return
    publication_id = publication_id.strip()
    doiless_publication_id = publication_id.split("doi:", 1)[-1]
    if not doiless_publication_id:
        lint_ctx.error("Empty DOI citation")
        return
    url = f"{base_url}/{doiless_publication_id}"
    r = requests.get(url)
    if r.status_code == 200:
        if publication_id != doiless_publication_id:
            lint_ctx.error("%s is valid, but Galaxy expects DOI without 'doi:' prefix" % publication_id)
        else:
            lint_ctx.info("%s is a valid DOI" % publication_id)
    elif r.status_code == 404:
        lint_ctx.error("%s is not a valid DOI" % publication_id)
    else:
        lint_ctx.warn("dx.doi returned unexpected status code %d" % r.status_code)


def lint_xsd(lint_ctx, schema_path, path):
    """Lint XML at specified path with supplied schema."""
    name = lint_ctx.object_name or os.path.basename(path)
    validator = validation.get_validator(require=True)
    validation_result = validator.validate(schema_path, path)
    if not validation_result.passed:
        msg = "Invalid XML found in file: %s. Errors [%s]"
        msg = msg % (name, validation_result.output)
        lint_ctx.error(msg)
    else:
        lint_ctx.info("File validates against XML schema.")


def _validate_doi_url(url, lint_ctx):
    """Validate DOI URL by checking CrossRef API."""
    match = re.match("https?://doi.org/(.*)$", url)
    if match is None:
        return False

    doi = match.group(1)
    xref_url = f"https://api.crossref.org/works/{doi}"
    return _validate_http_url(xref_url, lint_ctx=lint_ctx)


def _validate_http_url(url, lint_ctx, user_agent=None):
    """Validate HTTP/HTTPS URL."""
    headers = {"User-Agent": user_agent, "Accept": "*/*"} if user_agent else None
    r = None
    try:
        r = requests.get(url, headers=headers, stream=True, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        next(r.iter_content(1000))
        return True
    except Exception as e:
        if r is not None and r.status_code == 429:
            # too many requests
            return True
        elif r is not None and r.status_code in [403, 503] and "cloudflare" in r.text:
            # CloudFlare protection block
            return True
        else:
            lint_ctx.error(f"Error '{e}' accessing {url}")
            return False


def _validate_other_url(url, lint_ctx):
    """Validate non-HTTP URLs."""
    try:
        with urlopen(url) as handle:
            handle.read(100)
        return True
    except Exception as e:
        lint_ctx.error(f"Error '{e}' accessing {url}")
        return False


def lint_urls(root, lint_ctx):
    """Find referenced URLs and verify they are valid."""
    urls, docs = find_urls_for_xml(root)

    # This is from Google Chome on macOS, current at time of writing:
    BROWSER_USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 14_7_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/133.0.0.0 Safari/537.36"

    def validate_url(url, lint_ctx, user_agent=None):
        is_valid = False
        if re.match("https?://doi.org/(.*)$", url):
            is_valid = _validate_doi_url(url, lint_ctx)
        elif url.startswith("http://") or url.startswith("https://"):
            is_valid = _validate_http_url(url, lint_ctx, user_agent)
        else:
            is_valid = _validate_other_url(url, lint_ctx)

        if is_valid:
            lint_ctx.info("URL OK %s" % url)

    for url in urls:
        validate_url(url, lint_ctx)
    for url in docs:
        validate_url(url, lint_ctx, BROWSER_USER_AGENT)


__all__ = (
    "build_lint_args",
    "handle_lint_complete",
    "lint_dois",
    "lint_urls",
    "lint_xsd",
)
