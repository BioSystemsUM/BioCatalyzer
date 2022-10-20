import os
import shutil
from unittest import TestCase

import pandas as pd

from biocatalyzer.matcher import MSDataMatcher

from tests import TESTS_DATA_PATH


class MSDataMatcherTestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        self.new_output_folder = 'new_output_path/'
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

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

    def test_ms_data_matcher_properties_and_setters(self):
        ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data.tsv')
        compounds_to_match = os.path.join(TESTS_DATA_PATH, 'new_compounds_sample/new_compounds.tsv')
        ms = MSDataMatcher(ms_data_path=ms_data_path,
                           compounds_to_match=compounds_to_match,
                           output_path=self.output_folder,
                           mode='mass_diff',
                           tolerance=0.0015)

        output_path = ms.output_path
        self.assertEqual(output_path, self.output_folder)

        ms.output_path = self.new_output_folder
        shutil.rmtree(self.new_output_folder)

        with self.assertRaises(FileExistsError):
            ms.output_path = os.path.join(TESTS_DATA_PATH, 'results_sample/')

        ms.generate_ms_results()

        _ = ms.ms_data_path
        with self.assertRaises(FileNotFoundError):
            ms.ms_data_path = 'not_existing_path.tsv'

        ms.ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data_subsample.tsv')

        _ = ms.compounds_to_match
        with self.assertRaises(FileNotFoundError):
            ms.compounds_to_match = 'not_existing_path.tsv'

        ms.compounds_to_match = os.path.join(TESTS_DATA_PATH, 'new_compounds_sample/new_compounds_subsample.tsv')

        _ = ms.mode
        with self.assertRaises(ValueError):
            ms.mode = 'not_existing_mode'

        ms.mode = 'mass_diff'

        tl = ms.tolerance
        ms.tolerance = 0.0015 + tl

        _ = ms.matches
        with self.assertRaises(AttributeError):
            ms.matches = pd.DataFrame()

