#!/usr/bin/env python
import os
import sys
from xml.sax.saxutils import escape


def make_table(directory):
    ret = ['<table class="fileList">\n']
    for root, dirs, files in os.walk(directory):
        root = root.rstrip("/")
        ret.append(
            '   <tr><td class="directory">%s</td></tr>\n'
            % escape(os.path.split(root)[-1])
        )
        for file in files:
            ret.append(
                '   <tr><td class="file"><a href="%s">%s</a></td></tr>\n'
                % (os.path.join(os.path.split(root)[-1], file), escape(file))
            )
    ret.append("</table>")
    return "".join(ret)


def make_html(directory):
    return "\n".join(
        [
            "<html>" "<head>",
            "   <title>Search results</title>",
            '   <style type="text/css">',
            "      table.fileList { text-align: left; }",
            "      td.directory { font-weight: bold; }",
            "      td.file { padding-left: 4em; }",
            "   </style>",
            "</head>",
            "<body>",
            "<h1>Search Results</h1>",
            make_table(directory),
            "</body>",
            "</html>",
        ]
    )


if __name__ == "__main__":
    if len(sys.argv) == 2:
        directory_path = sys.argv[1]
    else:
        top = "."
    print(make_html(directory_path))
