# unit tests for general assets

import pytest

import blueswede


class TestVersion:
    def test_version_present(self):
        # any string len>0 returns true
        assert blueswede.__version__
