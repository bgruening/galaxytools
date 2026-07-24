import sys
from io import StringIO


class FieldStorage:
    def __init__(self, *args, **kwargs):
        self.value = self.file = self.name = self.filename = None
    def getvalue(self, key, default=None):
        return default
    def __bool__(self):
        return False


def escape(s, quote=True):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    if quote:
        s = s.replace('"', "&quot;").replace("'", "&#x27;")
    return s


def parse_header(line):
    return line, {}


def test():
    return StringIO()
