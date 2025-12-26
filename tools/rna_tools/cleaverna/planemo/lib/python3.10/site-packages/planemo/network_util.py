import socket
from http.client import BadStatusLine
from time import time as now
from urllib.error import URLError
from urllib.request import urlopen


def get_free_port() -> int:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("localhost", 0))
    port = sock.getsockname()[1]
    sock.close()
    return port


def wait_http_service(url, timeout=None):
    if timeout:
        end = now() + timeout

    while True:
        try:
            if timeout:
                next_timeout = end - now()
                if next_timeout < 0:
                    return False

            kwds = {} if timeout is None else dict(timeout=next_timeout)
            with urlopen(url, **kwds) as r:
                if r.getcode() != 200:
                    continue
            return True
        except OSError:
            pass
        except BadStatusLine:
            pass
        except URLError:
            pass


# code.activestate.com/recipes/576655-wait-for-network-service-to-appear
def wait_net_service(server, port, timeout=None):
    """Wait for network service to appear.

    :param int timeout: in seconds, if None or 0 wait forever
    :return: A ``bool`` - if ``timeout`` is ``None`` may return only ``True`` or
             throw an unhandled network exception.
    """
    if port is None:
        raise TypeError("wait_net_service passed NoneType port value.")

    port = int(port)

    if timeout:
        end = now() + timeout

    while True:
        s = socket.socket()
        # Following line prevents this method from interfering with process
        # it is waiting for on localhost.
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            if timeout:
                next_timeout = end - now()
                if next_timeout < 0:
                    return False
                else:
                    s.settimeout(next_timeout)

            s.connect((server, port))

        except socket.timeout:
            # this exception occurs only if timeout is set
            if timeout:
                return False

        except OSError:
            # if getattr(e, "errno") == 61:
            #    refused_connections += 1
            s.close()
        else:
            s.close()
            return True
