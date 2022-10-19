from unittest import TestCase

from biocatalyzer._utils import match_value


class TestUtils(TestCase):

    def test_match_value(self):
        self.assertTrue(match_value(10, [10.1, 10.2], 0.1)[0])
        self.assertFalse(match_value(10, [10.1, 10.2], 0.01)[0])
