import os
import re
from collections.abc import MutableMapping, MutableSequence
from typing import Any, AnyStr, Callable, Final, Optional, Union

import ruamel.yaml
from ruamel.yaml.comments import CommentedBase, CommentedMap, CommentedSeq

lineno_re: Final = re.compile(r"^(.*?:[0-9]+:[0-9]+: )(( *)(.*))")


def _add_lc_filename(r: ruamel.yaml.comments.CommentedBase, source: AnyStr) -> None:
    if isinstance(r, ruamel.yaml.comments.CommentedBase):
        r.lc.filename = source
    if isinstance(r, MutableSequence):
        for d in r:
            _add_lc_filename(d, source)
    elif isinstance(r, MutableMapping):
        for d in r.values():
            _add_lc_filename(d, source)


def relname(source: str) -> str:
    if source.startswith("file://"):
        source = source[7:]
        source = os.path.relpath(source)
    return source


def add_lc_filename(r: ruamel.yaml.comments.CommentedBase, source: str) -> None:
    _add_lc_filename(r, relname(source))


def reflow_all(text: str, maxline: Optional[int] = None) -> str:
    """Reflow text respecting a common prefix with line & col info."""
    if maxline is None:
        maxline = int(os.environ.get("COLUMNS", "100"))
    maxno = 0
    for line in text.splitlines():
        g = lineno_re.match(line)
        if not g:
            continue
        group = g.group(1)
        assert group is not None  # nosec
        maxno = max(maxno, len(group))
    maxno_text: Final = maxline - maxno
    msg: Final[list[str]] = []
    for line in text.splitlines():
        g = lineno_re.match(line)
        if not g:
            msg.append(line)
            continue
        pre = g.group(1)
        assert pre is not None  # nosec
        group2 = g.group(2)
        assert group2 is not None  # nosec
        reflowed = reflow(group2, maxno_text, g.group(3)).splitlines()
        msg.extend([pre.ljust(maxno, " ") + r for r in reflowed])
    return "\n".join(msg)


def reflow(text: str, maxline: int, shift: Optional[str] = "") -> str:
    """Reflow a single line of text."""
    maxline2: Final = max(maxline, 20)
    if len(text) > maxline2:
        sp = text.rfind(" ", 0, maxline2)
        if sp < 1:
            sp = text.find(" ", sp + 1)
            if sp == -1:
                sp = len(text)
        if sp < len(text):
            return f"{text[0:sp]}\n{shift}{reflow(text[sp + 1 :], maxline2, shift)}"
    return text


def _lineno(i: int, line: str, shift: str, bullet: str) -> str:
    if (r := lineno_re.match(line)) is not None:
        group1: Final = r.group(1)
        group2: Final = r.group(2)
        assert group1 is not None  # nosec
        assert group2 is not None  # nosec
        return group1 + (bullet if i == 0 else shift) + group2
    return (bullet if i == 0 else shift) + line


def indent(v: str, nolead: bool = False, shift: str = "  ", bullet: str = "  ") -> str:
    if nolead:
        return v.splitlines()[0] + "\n".join([shift + line for line in v.splitlines()[1:]])

    return "\n".join([_lineno(i, line, shift, bullet) for i, line in enumerate(v.splitlines())])


def bullets(textlist: list[str], bul: str) -> str:
    """Indent using the specified number of bullets."""
    if len(textlist) == 1:
        return textlist[0]
    return "\n".join(indent(t, bullet=bul) for t in textlist)


def strip_duplicated_lineno(text: str) -> str:
    """
    Strip duplicated line numbers.

    Same as :py:meth:`strip_dup_lineno` but without reflow.
    """
    pre: Optional[str] = None
    msg: Final = []
    for line in text.splitlines():
        g = lineno_re.match(line)
        if not g:
            msg.append(line)
            continue
        if g.group(1) != pre:
            msg.append(line)
            pre = g.group(1)
        else:
            group1 = g.group(1)
            group2 = g.group(2)
            assert group1 is not None  # nosec
            assert group2 is not None  # nosec
            msg.append(" " * len(group1) + group2)
    return "\n".join(msg)


def strip_dup_lineno(text: str, maxline: Optional[int] = None) -> str:
    """Strip duplicated line numbers."""
    if maxline is None:
        maxline = int(os.environ.get("COLUMNS", "100"))
    pre: Optional[str] = None
    msg: Final = []
    maxno = 0
    for line in text.splitlines():
        g = lineno_re.match(line)
        if not g:
            continue
        group1 = g.group(1)
        assert group1 is not None  # nosec
        maxno = max(maxno, len(group1))

    for line in text.splitlines():
        g = lineno_re.match(line)
        if not g:
            msg.append(line)
            continue
        if g.group(1) != pre:
            group3 = g.group(3)
            assert group3 is not None  # nosec
            shift = maxno + len(group3)
            group2 = g.group(2)
            assert group2 is not None  # nosec
            g2 = reflow(group2, maxline - shift, " " * shift)
            pre = g.group(1)
            assert pre is not None  # nosec
            msg.append(pre + " " * (maxno - len(pre)) + g2)
        else:
            group2 = g.group(2)
            assert group2 is not None  # nosec
            group3 = g.group(3)
            assert group3 is not None  # nosec
            g2 = reflow(group2, maxline - maxno, " " * (maxno + len(group3)))
            msg.append(" " * maxno + g2)
    return "\n".join(msg)


def cmap(
    d: Union[int, float, str, MutableMapping[str, Any], MutableSequence[Any], None],
    lc: Optional[list[int]] = None,
    fn: Optional[str] = None,
) -> Union[int, float, str, CommentedMap, CommentedSeq, None]:
    """
    Apply line+column & filename data through to the provided data.

    :param d: Target datastructure: number and strings will be passed through
              unchanged. Non-Commented container types with be transformed into
              the relevant Commented container types.
    :param lc: Line & Column information to be applied.
    :param fn: The Filename to store.
    :returns: The (transformed) datastructure.
    """
    if lc is None:
        lc = [0, 0, 0, 0]
    if fn is None:
        fn = "test"

    if isinstance(d, CommentedMap):
        fn = d.lc.filename if hasattr(d.lc, "filename") else fn
        for k, v in d.items():
            if d.lc.data is not None and k in d.lc.data:
                d[k] = cmap(v, lc=d.lc.data[k], fn=fn)
            else:
                d[k] = cmap(v, lc, fn=fn)
        return d
    if isinstance(d, CommentedSeq):
        fn = d.lc.filename if hasattr(d.lc, "filename") else fn
        for k2, v2 in enumerate(d):
            if d.lc.data is not None and k2 in d.lc.data:
                d[k2] = cmap(v2, lc=d.lc.data[k2], fn=fn)
            else:
                d[k2] = cmap(v2, lc, fn=fn)
        return d
    if isinstance(d, MutableMapping):
        cm: Final = CommentedMap()
        for k in sorted(d.keys()):
            v = d[k]
            if isinstance(v, CommentedBase):
                uselc = [v.lc.line, v.lc.col, v.lc.line, v.lc.col]
                vfn = v.lc.filename if hasattr(v.lc, "filename") else fn
            else:
                uselc = lc
                vfn = fn
            cm[k] = cmap(v, lc=uselc, fn=vfn)
            cm.lc.add_kv_line_col(k, uselc)
            cm.lc.filename = fn
        return cm
    if isinstance(d, MutableSequence):
        cs: Final = CommentedSeq()
        for k3, v3 in enumerate(d):
            if isinstance(v3, CommentedBase):
                uselc = [v3.lc.line, v3.lc.col, v3.lc.line, v3.lc.col]
                vfn = v3.lc.filename if hasattr(v3.lc, "filename") else fn
            else:
                uselc = lc
                vfn = fn
            cs.append(cmap(v3, lc=uselc, fn=vfn))
            cs.lc.add_kv_line_col(k3, uselc)
            cs.lc.filename = fn
        return cs
    return d


class SourceLine:
    def __init__(
        self,
        item: Any,
        key: Optional[Any] = None,
        raise_type: Callable[[str], Any] = str,
        include_traceback: bool = False,
    ) -> None:
        self.item = item
        self.key = key
        self.raise_type = raise_type
        self.include_traceback = include_traceback

    def __enter__(self) -> "SourceLine":
        return self

    def __exit__(
        self,
        exc_type: Any,
        exc_value: Any,
        tb: Any,
    ) -> None:
        if not exc_value:
            return
        raise self.makeError(str(exc_value)) from exc_value

    def file(self) -> Optional[str]:
        """Return the embedded filename."""
        if hasattr(self.item, "lc") and hasattr(self.item.lc, "filename"):
            return str(self.item.lc.filename)
        return None

    def start(self) -> Optional[tuple[int, int]]:
        """Determine the starting location."""
        if self.file() is None:
            return None
        if self.key is None or self.item.lc.data is None or self.key not in self.item.lc.data:
            return ((self.item.lc.line or 0) + 1, (self.item.lc.col or 0) + 1)
        return (
            (self.item.lc.data[self.key][0] or 0) + 1,
            (self.item.lc.data[self.key][1] or 0) + 1,
        )

    def end(self) -> Optional[tuple[int, int]]:
        """Empty, for now."""
        return None

    def makeLead(self) -> str:
        if self.file():
            lcol: Final = self.start()
            line, col = lcol if lcol else ("", "")
            return f"{self.file()}:{line}:{col}:"
        return ""

    def makeError(self, msg: str) -> Any:
        if not isinstance(self.item, ruamel.yaml.comments.CommentedBase):
            return self.raise_type(msg)
        errs: Final = []
        lead: Final = self.makeLead()
        for m in msg.splitlines():
            if bool(lineno_re.match(m)):
                errs.append(m)
            else:
                errs.append(f"{lead} {m}")
        return self.raise_type("\n".join(errs))
