"""
Tests on the GalaxyInstance object itself.
"""

import os
import time
import unittest

import pytest

from bioblend import ConnectionError
from bioblend.galaxy import GalaxyInstance
from . import test_util


class TestGalaxyInstance(unittest.TestCase):
    def setUp(self):
        # "connect" to a fake Galaxy instance
        self.gi = GalaxyInstance("http://localhost:56789", key="whatever")

    def test_url_attribute(self):
        assert self.gi.base_url == "http://localhost:56789"
        assert self.gi.url == "http://localhost:56789/api"
        # Test instance served at a subdirectory
        gi = GalaxyInstance("http://localhost:56789/galaxy/", key="whatever")
        assert gi.base_url == "http://localhost:56789/galaxy"
        assert gi.url == "http://localhost:56789/galaxy/api"

    def test_set_max_get_attempts(self):
        self.gi.max_get_attempts = 3
        assert 3 == self.gi.max_get_attempts

    def test_set_retry_delay(self):
        self.gi.get_retry_delay = 5.0
        assert 5.0 == self.gi.get_retry_delay

    def test_get_retry(self):
        # We set the client to try twice, with a delay of 5 seconds between
        # attempts. So, we expect the call to take at least 5 seconds before
        # failing.
        self.gi.max_get_attempts = 3
        self.gi.get_retry_delay = 2
        start = time.time()
        with pytest.raises(ConnectionError):
            self.gi.libraries.get_libraries()
        end = time.time()
        duration = end - start
        assert duration > self.gi.get_retry_delay * (self.gi.max_get_attempts - 1), "Didn't seem to retry long enough"

    def test_missing_scheme_fake_url(self):
        with pytest.raises(ValueError):
            GalaxyInstance("localhost:56789", key="whatever")

    @test_util.skip_unless_galaxy()
    def test_missing_scheme_real_url(self):
        galaxy_url = os.environ["BIOBLEND_GALAXY_URL"]
        # Strip the scheme from galaxy_url
        scheme_sep = "://"
        if scheme_sep in galaxy_url:
            galaxy_url = galaxy_url.partition(scheme_sep)[2]
        GalaxyInstance(url=galaxy_url)
