#!/usr/bin/env python3
# spec: simplest python web server with range support and multithreading that takes root path,
# port and bind address as command line arguments; by default uses the current dir as webroot,
# port 8000 and bind address of 0.0.0.0
# borrowed from https://github.com/danvk/RangeHTTPServer
# and reborrowed from https://gist.github.com/glowinthedark/b99900abe935e4ab4857314d647a9068
#
# The Apache 2.0 license copy in this repository is distributed with this code in accordance with that licence.
# https://www.apache.org/licenses/LICENSE-2.0.txt
# This part is not MIT licenced like the other components.

# APPENDIX: How to apply the Apache License to your work.

# To apply the Apache License to your work, attach the following
# boilerplate notice, with the fields enclosed by brackets "[]"
# replaced with your own identifying information. (Don't include
# the brackets!)  The text should be enclosed in the appropriate
# comment syntax for the file format. We also recommend that a
# file or class name and description of purpose be included on the
# same "printed page" as the copyright notice for easier
# identification within third-party archives.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse
import functools
import os
import re
import socketserver
import webbrowser
from http.server import SimpleHTTPRequestHandler


DEFAULT_PORT = 8081


def copy_byte_range(infile, outfile, start=None, stop=None, bufsize=16 * 1024):
    """Like shutil.copyfileobj, but only copy a range of the streams.

    Both start and stop are inclusive.
    """
    if start is not None:
        infile.seek(start)
    while 1:
        to_read = min(bufsize, stop + 1 - infile.tell() if stop else bufsize)
        buf = infile.read(to_read)
        if not buf:
            break
        outfile.write(buf)


BYTE_RANGE_RE = re.compile(r"bytes=(\d+)-(\d+)?$")


def parse_byte_range(byte_range):
    """Returns the two numbers in 'bytes=123-456' or throws ValueError.

    The last number or both numbers may be None.
    """
    if byte_range.strip() == "":
        return None, None

    m = BYTE_RANGE_RE.match(byte_range)
    if not m:
        raise ValueError("Invalid byte range %s" % byte_range)

    first, last = [x and int(x) for x in m.groups()]
    if last and last < first:
        raise ValueError("Invalid byte range %s" % byte_range)
    return first, last


class RangeRequestHandler(SimpleHTTPRequestHandler):
    """Adds support for HTTP 'Range' requests to SimpleHTTPRequestHandler

    The approach is to:
    - Override send_head to look for 'Range' and respond appropriately.
    - Override copyfile to only transmit a range when requested.
    """

    def handle(self):
        try:
            SimpleHTTPRequestHandler.handle(self)
        except Exception:
            # ignored, thrown whenever the client aborts streaming (broken pipe)
            pass

    def send_head(self):
        if "Range" not in self.headers:
            self.range = None
            return SimpleHTTPRequestHandler.send_head(self)
        try:
            self.range = parse_byte_range(self.headers["Range"])
        except ValueError:
            self.send_error(400, "Invalid byte range")
            return None
        first, last = self.range

        # Mirroring SimpleHTTPServer.py here
        path = self.translate_path(self.path)
        f = None
        ctype = self.guess_type(path)
        try:
            f = open(path, "rb")
        except IOError:
            self.send_error(404, "File not found")
            return None

        fs = os.fstat(f.fileno())
        file_len = fs[6]
        if first >= file_len:
            self.send_error(416, "Requested Range Not Satisfiable")
            return None

        self.send_response(206)
        self.send_header("Content-type", ctype)

        if last is None or last >= file_len:
            last = file_len - 1
        response_length = last - first + 1

        self.send_header("Content-Range", "bytes %s-%s/%s" % (first, last, file_len))
        self.send_header("Content-Length", str(response_length))
        self.send_header("Last-Modified", self.date_time_string(fs.st_mtime))
        self.end_headers()
        return f

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        return SimpleHTTPRequestHandler.end_headers(self)

    def copyfile(self, source, outputfile):
        if not self.range:
            return SimpleHTTPRequestHandler.copyfile(self, source, outputfile)

        # SimpleHTTPRequestHandler uses shutil.copyfileobj, which doesn't let
        # you stop the copying before the end of the file.
        start, stop = self.range  # set in send_head()
        copy_byte_range(source, outputfile, start, stop)


class ThreadedTCPServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    allow_reuse_address = True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tiny Python Web Server supporting range requests, for local viewing of unzipped Galaxy JBrowse2 configurations"
    )
    parser.add_argument(
        "--root",
        default=os.getcwd(),
        help="Root path to serve files from (default: current working directory)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to listen on (default: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--bind",
        default="127.0.0.1",
        help="IP address to bind to (default: 127.0.0.1 - use 0.0.0.0 to allow access on your network)",
    )
    args = parser.parse_args()

    handler = functools.partial(RangeRequestHandler, directory=args.root)

    webbrowser.open(f"http://{args.bind}:{args.port}")

    with ThreadedTCPServer((args.bind, args.port), handler) as httpd:
        print(
            f"Serving HTTP on {args.bind} port {args.port} (http://{args.bind}:{args.port}/)"
        )
        httpd.serve_forever()
