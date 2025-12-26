"""
Test ``schema-salad-doc``.

(also known as ``schema-salad-tool --print-doc``)

For convenience, tests are checking exact strings. In the event of changes in
the "mistune" package, makedoc.py, or other changes, feel free to modify the test
strings as long as the new HTML renders the same way in typical browsers.

Likewise, if the schema-salad metaschema changes and it is missing one or more
of the features tested below, then please copy those old features to a new file
and update the affected tests to use those new file(s).
"""

import hashlib
import inspect
import json
import tempfile
from io import StringIO
from pathlib import Path
from typing import Optional

import pytest

from schema_salad.makedoc import makedoc

from .util import get_data


def test_schema_salad_inherit_docs() -> None:
    """Test schema-salad-doc when types inherit and override values from parent types."""
    schema_path = get_data("tests/inherited-attributes.yml")
    stdout = StringIO()
    makedoc(stdout, schema_path)

    # The parent ID documentation (i.e. Parent ID) must appear exactly once.
    assert 1 == stdout.getvalue().count("Parent ID")


def generate_doc(schema_data: Optional[str] = None) -> str:
    """Avoid error when calling fixture directly."""
    stdout = StringIO()
    if schema_data:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".yml") as tmp_file:
            tmp_file.write(schema_data)
            tmp_file.flush()
            tmp_file.seek(0)
            makedoc(stdout, tmp_file.name)
    else:
        schema_path = get_data("metaschema/metaschema.yml")
        makedoc(stdout, schema_path)
    return stdout.getvalue()


@pytest.fixture(scope="session", name="metaschema_doc")
def fixture_metaschema_doc() -> str:
    """Pytest Fixture of the rendered HTML for the metaschema schema."""
    return generate_doc()


def test_doc_fenced_code_contents_preserved() -> None:
    """
    Fenced code contents are not interpreted as Markdown definitions and converted into erroneous HTML.

    An example of problem case is when a definition looks like a Markdown list
    (e.g.: a YAML array). It must not be converted into HTML contents with list tags.
    However, special characters (e.g.: ``<``, ``>``) must still be escaped,
    otherwise they will not be correctly rendered within an HTML ``<pre><code>`` block.
    """
    data = inspect.cleandoc(
        """
        Option one generic example:
        ```
        some_cwl_field:
          - key_field: a_complex_type1
            field2: foo
            field3: bar
          - key_field: a_complex_type2
            field2: foo2
            field3: bar2
          - key_field: a_complex_type3
        ```

        Option one generic example:
        ```
        inputs:
          - id: workflow_input01
            type: string
          - id: workflow_input02
            type: File
            format: http://edamontology.org/format_2572
        ```

        Special Characters:
        ```
        data:
          - test: value 1 < 2 is true
          - test: value 2 > 1 is false
        ```
        """
    )
    schema_data = json.dumps(
        [
            {
                "name": "https://w3id.org/cwl/salad#test",
                "type": "documentation",
                "doc": data,
            }
        ]
    )
    schema_doc = generate_doc(schema_data)
    schema_doc = schema_doc.split("</head>")[-1]  # just for nicer/readable assert diffs
    assert (
        "<p>Option one generic example:</p>\n"
        "<pre><code>some_cwl_field:\n"
        "  - key_field: a_complex_type1\n"
        "    field2: foo\n"
        "    field3: bar\n"
        "  - key_field: a_complex_type2\n"
        "    field2: foo2\n"
        "    field3: bar2\n"
        "  - key_field: a_complex_type3\n"
        "</code></pre>\n"
    ) in schema_doc
    assert (
        "<p>Option one generic example:</p>\n"
        "<pre><code>inputs:\n"
        "  - id: workflow_input01\n"
        "    type: string\n"
        "  - id: workflow_input02\n"
        "    type: File\n"
        "    format: http://edamontology.org/format_2572\n"
        "</code></pre>\n"
    ) in schema_doc
    assert (
        "<p>Special Characters:</p>\n"
        "<pre><code>data:\n"
        "  - test: value 1 &lt; 2 is true\n"
        "  - test: value 2 &gt; 1 is false\n"
        "</code></pre>\n"
    ) in schema_doc


def test_doc_headings_target_anchor(metaschema_doc: str) -> None:
    """Doc headers must have an id and section link."""
    assert (
        '<h1 id="Abstract" class="section">Abstract '
        '<a href="#Abstract">&sect;</a></h1>' in metaschema_doc
    )


def test_doc_render_table_of_contents(metaschema_doc: str) -> None:
    """The special Table of Contents token must be replaced with a rendered table."""
    assert "!--ToC--" not in metaschema_doc


def test_plain_links_autolinked(metaschema_doc: str) -> None:
    """Plan links should be treated as if they were wrapped in angle brackets."""
    assert (
        "This document is the product of the "
        '<a href="https://groups.google.com/forum/#!forum/common-workflow-language">'
        "Common Workflow Language working\ngroup</a>" in metaschema_doc
    )


def test_embedded_html_unescaped() -> None:
    """Raw HTML shouldn't get escaped."""
    schema_path = get_data("tests/inherited-attributes.yml")
    stdout = StringIO()
    makedoc(stdout, schema_path)
    html = stdout.getvalue()

    assert '<table class="table">' in html
    assert "&lt;table class=&quot;table&quot;&gt;" not in html


def test_multiline_list_entries_word_spacing(metaschema_doc: str) -> None:
    """Hanging indents in Markdown lists don't lead to wordsmushing."""
    assert "as itis poorly documented" not in metaschema_doc
    assert "base URI for the document used toresolve relative" not in metaschema_doc
    assert "The keys ofthe object are namespace prefixes" not in metaschema_doc
    assert "This field may list URIreferences to documents in RDF-XML" not in metaschema_doc
    assert "defines valid fields thatmake up a record type" not in metaschema_doc
    assert "set of symbols that arevalid value" not in metaschema_doc


def test_multiline_list_entries_without_indention(metaschema_doc: str) -> None:
    """Hanging indents are not required in Markdown lists."""
    # See https://daringfireball.net/projects/markdown/syntax#list
    # and https://spec.commonmark.org/0.30/#example-290
    # Some newlines in markdown are replaced by spaces purposely
    # to avoid invalid Markdown to HTML conversion as preserve words separated.
    assert (
        "<li><p>At least one record definition object which defines valid fields that\n"
        "make up a record type.  Record field definitions include the valid types\n"
        "that may be assigned to each field and annotations to indicate fields\n"
        'that represent identifiers and links, described below in "Semantic\n'
        'Annotations".</p>\n'
        "</li>\n"
        "<li><p>Any number of enumerated type objects which define a set of finite "
        "set of symbols that are\n"
        "valid value of the type.</p>\n"
        "</li>\n"
        "<li><p>Any number of documentation objects which allow in-line "
        "documentation of the schema.</p>\n"
        "</li>" in metaschema_doc
    )
    assert (
        "<li>At least one record definition object which defines valid fields that</li>"
        not in metaschema_doc
    )


def test_detect_changes_in_html(metaschema_doc: str, tmp_path: Path) -> None:
    """Catch all for changes in HTML output, please adjust if the changes are innocent."""
    # If the hash changed because the metaschema itself changed (without changes
    # to makedoc.py or the version of mistune) then you can directly update the
    # hash value below.
    #
    # Otherwise, follow this procedure verify that the changed HTML rendering
    # is acceptable (or use make 'check-metaschema-diff' and 'compute-metaschema-hash'):
    # 1. Render the metaschema schema into using an older, known-working version
    #    of schema-salad:
    #    `schema-salad-doc schema_salad/metaschema/metaschema.yml > /tmp/metaschema.orig.html`
    # 2. Render the metaschema schema into using proposed changed codebase
    #    `schema-salad-doc schema_salad/metaschema/metaschema.yml > /tmp/metaschema.new.html`
    # 3. Confirm the other tests in this file pass using the new code/mistune,
    #    adjusting the test strings if the changes are truly innocent.
    # 4. Check the `diff` between the saved HTML pages to check for obvious problems
    #    `diff /tmp/metaschema.orig.html /tmp/metaschema.new.html`
    # 5. Check the HTML in both Firefox and Chrome, especially near areas
    #    of differences in the diff
    # 6. If the changes are agreeable, then update the hash below
    hasher = hashlib.sha256()
    hasher.update(metaschema_doc.encode("utf-8"))
    result = tmp_path / "result.html"
    with open(result, "w") as h:
        h.write(metaschema_doc)
    assert (
        hasher.hexdigest() == "bfab566e522ea2955bd358c071c8ea9b07591fec7949a71ae98ef3c7d11d5a7b"
    ), result
