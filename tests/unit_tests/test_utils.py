from unittest import TestCase

import numpy as np
import pandas as pd

from biocatalyzer._utils import match_value, _empty_dfs, _merge_fields


class TestUtils(TestCase):

    def test_match_value(self):
        self.assertEqual(len(match_value(10, [10.1, 10.2], 0.1)), 1)
        self.assertEqual(len(match_value(10.1, [10, 10.2], 0.1)), 2)
        self.assertEqual(len(match_value(10, [10.1, 10.2], 0.01)), 0)
        self.assertEqual(len(match_value(10, [10.1, 10.3], 0.1)), 1)

    def test_empty_dfs(self):
        dfs = [pd.DataFrame(), pd.DataFrame()]
        self.assertTrue(_empty_dfs(dfs))
        dfs = [pd.DataFrame(), pd.DataFrame({'a': [1]})]
        self.assertFalse(_empty_dfs(dfs))

    def test_merge_fields(self):
        self.assertEqual('1.1.1.1;1.1.1.2', _merge_fields('1.1.1.1;1.1.1.2;;1.1.1.1'))
        self.assertEqual('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4', _merge_fields('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4'))
        self.assertEqual('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4',
                         _merge_fields('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4;1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4'))
        self.assertTrue(np.isnan(_merge_fields(';;;')))

    def test_merge_fields2(self):
        self.assertEqual('field', _merge_fields('field'))
        self.assertEqual('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4', _merge_fields('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4'))
        self.assertEqual('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4',
                         _merge_fields('1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4;1.1.1.1;1.1.1.2;1.1.1.3;1.1.1.4'))
        self.assertTrue(np.isnan(_merge_fields('')))

