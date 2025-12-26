import argparse
import copy
import html
import logging
import os
import re
import sys
from collections.abc import MutableMapping, MutableSequence
from io import StringIO, TextIOWrapper
from typing import IO, Any, Final, Literal, Optional, Union, cast
from urllib.parse import urldefrag

from mistune import create_markdown
from mistune.markdown import Markdown
from mistune.renderers.html import HTMLRenderer

from .exceptions import SchemaSaladException, ValidationException
from .schema import avro_field_name, extend_and_specialize, get_metaschema
from .utils import add_dictlist, aslist
from .validate import avro_type_name

PluginName = Literal[
    "url",
    "strikethrough",
    "footnotes",
    "table",
    "task_lists",
    "def_list",
    "abbr",
]

_logger: Final = logging.getLogger("salad")

fenced_code_pattern: Final = re.compile(
    r"( *)(`{3,}|~{3,})([^`\n]*)\n(?:|([\s\S]*?)\n)(?:\1\2[~`]* *\n+|$)"
)
""" Pattern inspired from 'mistune.block_parser.BlockParser.FENCED_CODE'.
However, instead of the initial ' {0,3}' part to match any indented fenced-code,
use any quantity of spaces, as long as they match at the end as well (using '\1').
Because of nested fenced-code in lists, it can be more indented than "normal"."""


def escape_html(s: str) -> str:
    """Escape HTML but otherwise preserve single quotes."""
    return html.escape(html.unescape(s)).replace("&#x27;", "'")


def vocab_type_name(url: str) -> str:
    """Remove the avro namespace, if any."""
    return avro_type_name(url).split(".")[-1]


def has_types(items: Any) -> list[str]:
    """Retrieve all the types of a record."""
    r: Final[list[str]] = []
    if isinstance(items, MutableMapping):
        if items["type"] == "https://w3id.org/cwl/salad#record":
            return [items["name"]]
        for n in ("type", "items", "values"):
            if n in items:
                r.extend(has_types(items[n]))
        return r
    if isinstance(items, MutableSequence):
        for i in items:
            r.extend(has_types(i))
        return r
    if isinstance(items, str):
        return [items]
    return []


def linkto(item: str) -> str:
    frg: Final = urldefrag(item)[1]
    return f"[{frg}](#{to_id(frg)})"


class MyRenderer(HTMLRenderer):
    """Custom renderer with different representations of selected HTML tags."""

    def heading(self, text: str, level: int, **attrs: Any) -> str:
        """Override HTML heading creation with text IDs."""
        return """<h{} id="{}" class="section">{} <a href="#{}">&sect;</a></h{}>""".format(
            level, to_id(text), text, to_id(text), level
        )

    def text(self, text: str) -> str:
        """Don't escape quotation marks."""
        # avoid convert of & if already escaped
        text = html.unescape(text)
        # html.escape does both single/double quotes ('/")
        # mistune.util.escape does only double quotes
        return html.escape(text, quote=self._escape)

    def inline_html(self, html: str) -> str:
        """Don't escape characters in predefined HTML within paragraph tags."""
        return html + "\n"

    def block_html(self, html: str) -> str:
        """Don't escape characters nor wrap predefined HTML within paragraph tags."""
        return html + "\n"

    def block_code(self, code: str, info: Optional[str] = None) -> str:
        """Don't escape quotation marks."""
        text = "<pre><code"
        if info is not None:
            info = info.strip()
        if info:
            lang = info.split(None, 1)[0]
            lang = escape_html(lang)
            text += ' class="language-' + lang + '"'
        return text + ">" + html.escape(code, quote=self._escape) + "</code></pre>\n"


def patch_fenced_code(original_markdown_text: str, modified_markdown_text: str) -> str:
    """Reverts fenced code fragments found in the modified contents back to their original definition."""
    matches_original: Final = list(re.finditer(fenced_code_pattern, original_markdown_text))
    matches_modified: Final = list(re.finditer(fenced_code_pattern, modified_markdown_text))
    if len(matches_original) != len(matches_modified):
        raise ValueError("Cannot patch fenced code definitions with inconsistent matches.")
    result = ""
    begin = 0
    for original, modified in zip(matches_original, matches_modified):
        ori_s, ori_e = original.start(), original.end()
        mod_s, mod_e = modified.start(), modified.end()
        result += modified_markdown_text[begin:mod_s]  # add text in between matches
        result += original_markdown_text[ori_s:ori_e]  # revert the fenced code
        begin = mod_e  # skip over the modified fenced code for next match
    result += modified_markdown_text[begin:]  # left over text after last match
    return result


def to_id(text: str) -> str:
    textid = text
    if text[0] in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"):
        try:
            textid = text[text.index(" ") + 1 :]
        except ValueError:
            pass
    return textid.replace(" ", "_")


class ToC:
    def __init__(self) -> None:
        self.first_toc_entry = True
        self.numbering = [0]
        self.toc = ""
        self.start_numbering = True

    def add_entry(self, thisdepth: int, title: str) -> str:
        """Add an entry to the table of contents."""
        depth: Final = len(self.numbering)
        if thisdepth < depth:
            self.toc += "</ol>"
            for _ in range(0, depth - thisdepth):
                self.numbering.pop()
                self.toc += "</li></ol>"
            self.numbering[-1] += 1
        elif thisdepth == depth:
            if not self.first_toc_entry:
                self.toc += "</ol>"
            else:
                self.first_toc_entry = False
            self.numbering[-1] += 1
        elif thisdepth > depth:
            self.numbering.append(1)

        num: Final = (
            "{}.{}".format(self.numbering[0], ".".join([str(n) for n in self.numbering[1:]]))
            if self.start_numbering
            else ""
        )
        self.toc += f"""<li><a href="#{to_id(title)}">{num} {title}</a><ol>\n"""
        return num

    def contents(self, idn: str) -> str:
        toc = """<h1 id="{}">Table of contents</h1>
               <nav class="tocnav"><ol>{}""".format(
            idn, self.toc
        )
        toc += "</ol>"
        for _ in range(0, len(self.numbering)):
            toc += "</li></ol>"
        toc += """</nav>"""
        return toc


basicTypes: Final = (
    "https://w3id.org/cwl/salad#null",
    "http://www.w3.org/2001/XMLSchema#boolean",
    "http://www.w3.org/2001/XMLSchema#int",
    "http://www.w3.org/2001/XMLSchema#long",
    "http://www.w3.org/2001/XMLSchema#float",
    "http://www.w3.org/2001/XMLSchema#double",
    "http://www.w3.org/2001/XMLSchema#string",
    "https://w3id.org/cwl/salad#record",
    "https://w3id.org/cwl/salad#enum",
    "https://w3id.org/cwl/salad#array",
)


def number_headings(toc: ToC, maindoc: str) -> str:
    mdlines = []
    skip = False
    for line in maindoc.splitlines():
        if line.strip() == "# Introduction":
            toc.start_numbering = True
            toc.numbering.clear()
            toc.numbering.append(0)

        if "```" in line:
            skip = not skip

        if not skip:
            m = re.match(r"^(#+) (.*)", line)
            if m is not None:
                group1 = m.group(1)
                assert group1 is not None  # nosec
                group2 = m.group(2)
                assert group2 is not None  # nosec
                num = toc.add_entry(len(group1), group2)
                line = f"{group1} {num} {group2}"
            line = re.sub(r"^(https?://\S+)", r"[\1](\1)", line)
        mdlines.append(line)

    return "\n".join(mdlines)


def fix_doc(doc: Union[list[str], str]) -> str:
    """Concatenate doc strings, replacing email addresses with mailto links."""
    docstr: Final = "".join(doc) if isinstance(doc, MutableSequence) else doc
    return "\n".join(
        [re.sub(r"<([^>@]+@[^>]+)>", r"[\1](mailto:\1)", d) for d in docstr.splitlines()]
    )


def _extendsfrom(
    item: dict[str, Any], ex: list[dict[str, Any]], typemap: dict[str, dict[str, str]]
) -> None:
    if "extends" in item:
        for e in aslist(item["extends"]):
            ex.insert(0, typemap[e])
            _extendsfrom(typemap[e], ex, typemap)


class RenderType:
    def __init__(
        self,
        toc: ToC,
        j: list[dict[str, Any]],
        renderlist: list[str],
        redirects: dict[str, str],
        primitiveType: str,
    ) -> None:
        self.typedoc = StringIO()
        self.toc: Final = toc
        self.subs: Final[dict[str, str]] = {}
        self.docParent: Final[dict[str, list[str]]] = {}
        self.docAfter: Final[dict[str, list[str]]] = {}
        self.rendered: Final[set[str]] = set()
        self.redirects: Final = redirects
        self.title: Optional[str] = None
        self.primitiveType: Final = primitiveType

        for t in j:
            if "extends" in t:
                for e in aslist(t["extends"]):
                    add_dictlist(self.subs, e, t["name"])
                    # if "docParent" not in t and "docAfter" not in t:
                    #    add_dictlist(self.docParent, e, t["name"])

            if t.get("docParent"):
                add_dictlist(self.docParent, t["docParent"], t["name"])

            if t.get("docChild"):
                for c in aslist(t["docChild"]):
                    add_dictlist(self.docParent, t["name"], c)

            if t.get("docAfter"):
                add_dictlist(self.docAfter, t["docAfter"], t["name"])

        metaschema_loader: Final = get_metaschema()[2]
        alltypes: Final = extend_and_specialize(j, metaschema_loader)

        self.typemap: Final[dict[str, dict[str, str]]] = {}
        self.uses: Final[dict[str, list[tuple[str, str]]]] = {}
        self.record_refs: Final[dict[str, list[str]]] = {}
        for entry in alltypes:
            self.typemap[entry["name"]] = entry
            try:
                if entry["type"] == "record":
                    self.record_refs[entry["name"]] = []
                    fields: Union[str, list[dict[str, str]]] = entry.get("fields", [])
                    if isinstance(fields, str):
                        raise KeyError("record fields must be a list of mappings")
                    else:
                        for f in fields:
                            p = has_types(f)
                            for tp in p:
                                if tp not in self.uses:
                                    self.uses[tp] = []
                                if (entry["name"], f["name"]) not in self.uses[tp]:
                                    _, frg1 = urldefrag(t["name"])
                                    _, frg2 = urldefrag(f["name"])
                                    self.uses[tp].append((frg1, frg2))
                                if (
                                    tp not in basicTypes
                                    and tp not in self.record_refs[entry["name"]]
                                ):
                                    self.record_refs[entry["name"]].append(tp)
            except KeyError:
                _logger.error("Did not find 'type' in %s", t)
                _logger.error("record refs is %s", self.record_refs)
                raise

        for entry in alltypes:
            if entry["name"] in renderlist or (
                (not renderlist)
                and ("extends" not in entry)
                and ("docParent" not in entry)
                and ("docAfter" not in entry)
            ):
                self.render_type(entry, 1)

    def typefmt(
        self,
        tp: Any,
        redirects: dict[str, str],
        nbsp: bool = False,
        jsonldPredicate: Optional[Union[dict[str, str], str]] = None,
    ) -> str:
        if isinstance(tp, MutableSequence):
            if nbsp and len(tp) <= 3:
                return "&nbsp;|&nbsp;".join(
                    [self.typefmt(n, redirects, jsonldPredicate=jsonldPredicate) for n in tp]
                )
            return " | ".join(
                [self.typefmt(n, redirects, jsonldPredicate=jsonldPredicate) for n in tp]
            )
        if isinstance(tp, MutableMapping):
            if tp["type"] == "https://w3id.org/cwl/salad#array":
                ar = "array&lt;{}&gt;".format(self.typefmt(tp["items"], redirects, nbsp=True))
                if isinstance(jsonldPredicate, dict) and "mapSubject" in jsonldPredicate:
                    if "mapPredicate" in jsonldPredicate:
                        ar += " | "
                        if len(ar) > 40:
                            ar += "<br>"

                        ar += (
                            "<a href='#map'>map</a>&lt;<code>{}</code>"
                            ",&nbsp;<code>{}</code> | {}&gt".format(
                                jsonldPredicate["mapSubject"],
                                jsonldPredicate["mapPredicate"],
                                self.typefmt(tp["items"], redirects),
                            )
                        )
                    else:
                        ar += " | "
                        if len(ar) > 40:
                            ar += "<br>"
                        ar += "<a href='#map'>map</a>&lt;<code>{}</code>,&nbsp;{}&gt".format(
                            jsonldPredicate["mapSubject"],
                            self.typefmt(tp["items"], redirects),
                        )
                return ar
            if tp["type"] in (
                "https://w3id.org/cwl/salad#record",
                "https://w3id.org/cwl/salad#enum",
            ):
                frg: Final = vocab_type_name(tp["name"])
                if tp["name"] in redirects:
                    return """<a href="{}">{}</a>""".format(redirects[tp["name"]], frg)
                if tp["name"] in self.typemap:
                    return f"""<a href="#{to_id(frg)}">{frg}</a>"""
                if tp["type"] == "https://w3id.org/cwl/salad#enum" and len(tp["symbols"]) == 1:
                    return "constant value <code>{}</code>".format(
                        avro_field_name(tp["symbols"][0])
                    )
                return frg
            if isinstance(tp["type"], MutableMapping):
                return self.typefmt(tp["type"], redirects)
        else:
            if str(tp) in redirects:
                return f"""<a href="{redirects[tp]}">{redirects[tp]}</a>"""  # noqa: B907
            if str(tp) in basicTypes:
                return """<a href="{}">{}</a>""".format(
                    self.primitiveType, vocab_type_name(str(tp))
                )
            frg2: Final = urldefrag(tp)[1]
            if frg2 != "":
                tp = frg2
            return f"""<a href="#{to_id(tp)}">{tp}</a>"""
        raise SchemaSaladException("We should not be here!")

    def render_type(self, f: dict[str, Any], depth: int) -> None:
        """Render a type declaration."""
        if f["name"] in self.rendered or f["name"] in self.redirects:
            return
        self.rendered.add(f["name"])

        if f.get("abstract"):
            return

        if "doc" not in f:
            f["doc"] = ""

        f["type"] = copy.deepcopy(f)
        f["doc"] = ""
        f = f["type"]

        if "doc" not in f:
            f["doc"] = ""

        ex: Final = [f]
        _extendsfrom(f, ex, self.typemap)

        enumDesc = {}
        if f["type"] == "enum" and isinstance(f["doc"], MutableSequence):
            for e in ex:
                for i in e["doc"]:
                    idx = i.find(":")
                    if idx > -1:
                        enumDesc[i[:idx]] = i[idx + 1 :]
                e["doc"] = [i for i in e["doc"] if i.find(":") == -1 or i.find(" ") < i.find(":")]

        f["doc"] = fix_doc(f["doc"])

        if f["type"] == "record":
            for field in f.get("fields", []):
                if "doc" not in field:
                    field["doc"] = ""

        if f["type"] != "documentation":
            lines = []
            for line in f["doc"].splitlines():
                if len(line) > 0 and line[0] == "#":
                    line = ("#" * depth) + line
                lines.append(line)
            f["doc"] = "\n".join(lines)

            frg = urldefrag(f["name"])[1]
            num = self.toc.add_entry(depth, frg)
            doc = "{} {} {}\n".format(("#" * depth), num, frg)
        else:
            doc = ""

        # Save the first line of the first type definition for the
        # HTML <title> tag
        if self.title is None and f["doc"]:
            self.title = f["doc"].partition("\n")[0]
            if self.title.startswith("# "):
                self.title = self.title[2:]

        if f["type"] == "documentation":
            f["doc"] = number_headings(self.toc, f["doc"])

        doc = doc + "\n\n" + f["doc"]
        plugins: Final[list[PluginName]] = [
            "strikethrough",
            "footnotes",
            "table",
            "url",
        ]
        # if escape active, wraps literal HTML into '<p> {HTML} </p>'
        # we must pass it to both since 'MyRenderer' is predefined
        escape = False
        markdown2html: Markdown = create_markdown(
            renderer=MyRenderer(escape=escape),
            plugins=plugins,
            escape=escape,
        )
        doc = cast(str, markdown2html(doc))

        if f["type"] == "record":
            doc += "<h3>Fields</h3>"
            doc += """
<div class="responsive-table">
<div class="row responsive-table-header">
<div class="col-xs-3 col-lg-2">field</div>
<div class="col-xs-2 col-lg-1">required</div>
<div class="col-xs-7 col-lg-3">type</div>
<div class="col-xs-12 col-lg-6 description-header">description</div>
</div>"""
            required = []
            optional = []
            for i in f.get("fields", []):
                tp = i["type"]
                if isinstance(tp, MutableSequence) and tp[0] == "https://w3id.org/cwl/salad#null":
                    opt = False
                    tp = tp[1:]
                else:
                    opt = True

                desc = i["doc"]

                rfrg = avro_field_name(i["name"])
                tr = """
<div class="row responsive-table-row">
<div class="col-xs-3 col-lg-2"><code>{}</code></div>
<div class="col-xs-2 col-lg-1">{}</div>
<div class="col-xs-7 col-lg-3">{}</div>
<div class="col-xs-12 col-lg-6 description-col">{}</div>
</div>""".format(
                    rfrg,
                    "required" if opt else "optional",
                    self.typefmt(tp, self.redirects, jsonldPredicate=i.get("jsonldPredicate")),
                    markdown2html(desc),
                )
                if opt:
                    required.append(tr)
                else:
                    optional.append(tr)
            for i in required + optional:
                doc += i
            doc += """</div>"""
        elif f["type"] == "enum":
            doc += "<h3>Symbols</h3>"
            doc += """<table class="table table-striped">"""
            doc += "<tr><th>symbol</th><th>description</th></tr>"
            for e in ex:
                for i in e.get("symbols", []):
                    doc += "<tr>"
                    efrg = avro_field_name(i)
                    doc += "<td><code>{}</code></td><td>{}</td>".format(
                        efrg, enumDesc.get(efrg, "")
                    )
                    doc += "</tr>"
            doc += """</table>"""
        f["doc"] = doc

        self.typedoc.write(f["doc"])

        subs = self.docParent.get(f["name"], []) + self.record_refs.get(f["name"], [])
        if len(subs) == 1:
            self.render_type(self.typemap[subs[0]], depth)
        else:
            for s in subs:
                self.render_type(self.typemap[s], depth + 1)

        for s in self.docAfter.get(f["name"], []):
            self.render_type(self.typemap[s], depth)


def avrold_doc(
    j: list[dict[str, Any]],
    outdoc: IO[Any],
    renderlist: list[str],
    redirects: dict[str, str],
    brand: str,
    brandlink: str,
    primtype: str,
    brandstyle: Optional[str] = None,
    brandinverse: Optional[bool] = False,
) -> None:
    toc: Final = ToC()
    toc.start_numbering = False

    rt: Final = RenderType(toc, j, renderlist, redirects, primtype)
    content = rt.typedoc.getvalue()

    if brandstyle is None:
        bootstrap_url = "https://cdn.jsdelivr.net/npm/bootstrap@5.0.2/dist/css/bootstrap.min.css"
        bootstrap_integrity = (
            "sha384-EVSTQN3/azprG1Anm3QDgpJLIm9Nao0Yz1ztcQTwFspd3yD65VohhpuuCOmLASjC"
        )
        brandstyle_template = '<link rel="stylesheet" href={} integrity={} crossorigin="anonymous">'
        brandstyle = brandstyle_template.format(bootstrap_url, bootstrap_integrity)

    picturefill_url: Final = (
        "https://cdn.rawgit.com/scottjehl/picturefill/3.0.2/dist/picturefill.min.js"
    )
    picturefill_integrity: Final = (
        "sha384-ZJsVW8YHHxQHJ+SJDncpN90d0EfAhPP+yA94n+EhSRzhcxfo84yMnNk+v37RGlWR"
    )
    outdoc.write(
        """
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    {}
    <script>
    // Picture element HTML5 shiv
    document.createElement( "picture" );
    </script>
    <script src="{}"
        integrity="{}"
        crossorigin="anonymous" async></script>
    """.format(
            brandstyle, picturefill_url, picturefill_integrity
        )
    )

    outdoc.write(f"<title>{rt.title}</title>")

    outdoc.write(
        """
    <style>
    :target {
      padding-top: 61px;
      margin-top: -61px;
    }
    .tocnav ol {
      list-style: none
    }
    pre {
      margin: 0 2em 10px 2em;
      padding: 9.5px;
      line-height: 1.42857143;
      color: #333;
      word-break: break-all;
      word-wrap: break-word;
      background-color: #f5f5f5;
      border: 1px solid #ccc;
      border-radius: 4px;
    }
    pre code {
      padding: 0;
      font-size: inherit;
      color: inherit;
      white-space: pre-wrap;
      background-color: transparent;
      border-radius: 0;
    }
    code {
      background-color: #f9f2f4;
      border-radius: 4px;
      padding: 2px 4px;
      color: #c7254e;
    }
    blockquote {
      padding: 10px 20px 1px 20px;
      border-left: 5px solid #eee;
    }
    a {
      text-decoration: none;
    }
    a code {
      color: #c7254e;
    }
    .section a {
      visibility: hidden;
    }
    .section:hover a {
      visibility: visible;
      color: rgb(201, 201, 201);
    }
    .responsive-table-header {
      text-align: left;
      padding: 8px;
      vertical-align: top;
      font-weight: bold;
      border-top-color: rgb(221, 221, 221);
      border-top-style: solid;
      border-top-width: 1px;
      background-color: #f9f9f9
    }
    .responsive-table > .responsive-table-row {
      text-align: left;
      padding: 8px;
      vertical-align: top;
      border-top-color: rgb(221, 221, 221);
      border-top-style: solid;
      border-top-width: 1px;
    }
    @media (min-width: 0px), print {
      .description-header {
        display: none;
      }
      .description-col {
        margin-top: 1em;
        margin-left: 1.5em;
      }
    }
    @media (min-width: 1170px) {
      .description-header {
        display: inline;
      }
      .description-col {
        margin-top: 0px;
        margin-left: 0px;
      }
    }
    .responsive-table-row:nth-of-type(odd) {
       background-color: #f9f9f9
    }
    </style>
    </head>
    <body>
    """
    )

    navbar_extraclass: Final = "navbar-inverse" if brandinverse else ""
    outdoc.write(
        """
      <nav class="navbar sticky-top navbar-expand-lg navbar-light bg-light {}">
        <div class="container">
          <a class="navbar-brand" href="{}">{}</a>
    """.format(
            navbar_extraclass, brandlink, brand
        )
    )

    if "<!--ToC-->" in content:
        content = content.replace("<!--ToC-->", toc.contents("toc"))
        outdoc.write(
            """
              <ul class="navbar-nav me-auto">
                <li class="nav-item"><a class="nav-link" href="#toc">Table of contents</a></li>
              </ul>
        """
        )

    outdoc.write(
        """
        </div>
      </nav>
    """
    )

    outdoc.write(
        """
    <div class="container mt-4">
    """
    )

    outdoc.write(
        """
    <div class="row">
    """
    )

    outdoc.write(
        """
    <div class="col-md-12" role="main" id="main">"""
    )

    outdoc.write(content)

    outdoc.write("""</div>""")

    outdoc.write(
        """
    </div>
    </div>
    </body>
    </html>"""
    )


def arg_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser: Final = argparse.ArgumentParser()
    parser.add_argument("schema")
    parser.add_argument("--only", action="append")
    parser.add_argument("--redirect", action="append")
    parser.add_argument("--brand")
    parser.add_argument("--brandlink")
    parser.add_argument("--brandstyle")
    parser.add_argument("--brandinverse", default=False, action="store_true")
    parser.add_argument("--primtype", default="#PrimitiveType")
    parser.add_argument("--debug", action="store_true")
    return parser


def main() -> None:
    """Shortcut entrypoint."""
    args: Final = arg_parser().parse_args()
    if args.debug:
        _logger.setLevel(logging.DEBUG)
    makedoc(
        sys.stdout,
        args.schema,
        args.redirect,
        args.only,
        args.brand,
        args.brandlink,
        args.primtype,
        args.brandstyle,
        args.brandinverse,
    )


def makedoc(
    stdout: IO[Any],
    schema: str,
    redirects: Optional[list[str]] = None,
    only: Optional[list[str]] = None,
    brand: Optional[str] = None,
    brandlink: Optional[str] = None,
    primtype: Optional[str] = None,
    brandstyle: Optional[str] = None,
    brandinverse: Optional[bool] = False,
) -> None:
    """Emit HTML representation of a given schema."""
    s: Final[list[dict[str, Any]]] = []
    with open(schema, encoding="utf-8") as f:
        if schema.endswith("md"):
            s.append(
                {
                    "name": os.path.splitext(os.path.basename(schema))[0],
                    "type": "documentation",
                    "doc": f.read(),
                }
            )
        else:
            uri = "file://" + os.path.abspath(schema)
            metaschema_loader = get_metaschema()[2]
            j = metaschema_loader.resolve_ref(uri, "")[0]
            if isinstance(j, MutableSequence):
                s.extend(j)
            elif isinstance(j, MutableMapping):
                s.append(j)
            else:
                raise ValidationException("Schema must resolve to a list or a dict")
    redirect = {}
    for r in redirects or []:
        redirect[r.split("=")[0]] = r.split("=")[1]
    renderlist = only or []
    if hasattr(stdout, "buffer") and getattr(stdout, "encoding", None) != "UTF-8":
        wrapped_stdout: IO[Any] = TextIOWrapper(stdout.buffer, encoding="utf-8")
    else:
        wrapped_stdout = stdout
    avrold_doc(
        s,
        wrapped_stdout,
        renderlist,
        redirect,
        brand or "",
        brandlink or "",
        primtype or "",
        brandstyle=brandstyle,
        brandinverse=brandinverse,
    )


if __name__ == "__main__":
    main()
