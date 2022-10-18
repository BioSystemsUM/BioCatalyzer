import os
import shutil
from unittest import TestCase

import pandas as pd

from biocatalyzer.matcher import MSDataMatcher

from tests import TESTS_DATA_PATH


class MSDataMatcherTestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        os.mkdir(self.output_folder)

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestMSDataMatcher(MSDataMatcherTestCase, TestCase):

    def test_ms_data_matcher_mass_mode(self):
        ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data.tsv')
        compounds_to_match = os.path.join(TESTS_DATA_PATH, 'new_compounds_sample/new_compounds.tsv')
        ms = MSDataMatcher(ms_data_path=ms_data_path,
                           compounds_to_match=compounds_to_match,
                           output_path=self.output_folder,
                           mode='mass',
                           tolerance=0.0015)

        ms.generate_ms_results()

        self.assertEqual(ms.mode, 'mass')
        self.assertEqual(ms.tolerance, 0.0015)
        self.assertEqual(ms.compounds_to_match.shape, (266, 8))
        self.assertIsInstance(ms.matches, pd.DataFrame)
        self.assertEqual(ms.matches.shape, (0, 8))

    def test_ms_data_matcher_massdiff_mode(self):
        ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data.tsv')
        compounds_to_match = os.path.join(TESTS_DATA_PATH, 'new_compounds_sample/new_compounds.tsv')
        ms = MSDataMatcher(ms_data_path=ms_data_path,
                           compounds_to_match=compounds_to_match,
                           output_path=self.output_folder,
                           mode='mass_diff',
                           tolerance=0.0015)

        ms.generate_ms_results()

        self.assertEqual(ms.mode, 'mass_diff')
        self.assertEqual(ms.tolerance, 0.0015)
        self.assertEqual(ms.compounds_to_match.shape, (266, 9))
        self.assertIsInstance(ms.matches, pd.DataFrame)
        self.assertEqual(ms.matches.shape, (0, 9))

