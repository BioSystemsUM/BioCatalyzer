import shutil
from pathlib import Path
from unittest import TestCase

import pandas as pd

from biocatalyzer.matcher import MSDataMatcher

from tests import TESTS_DATA_PATH


class MSDataMatcherTestCase(TestCase):

    def setUp(self):
        self.output_folder = TESTS_DATA_PATH / 'results_sample'
        self.new_output_folder = TESTS_DATA_PATH / 'new_results_sample'
        # Ensure the directories exist
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.new_output_folder.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        if self.output_folder.exists():
            shutil.rmtree(self.output_folder)
        if self.new_output_folder.exists():
            shutil.rmtree(self.new_output_folder)


class TestMSDataMatcher(MSDataMatcherTestCase, TestCase):

    def test_ms_data_matcher(self):
        ms_data_path = TESTS_DATA_PATH / 'ms_data_sample' / 'ms_data.tsv'
        compounds_to_match = TESTS_DATA_PATH / 'new_compounds_sample' / 'new_compounds.tsv'
        ms = MSDataMatcher(ms_data_path=ms_data_path.as_posix(),
                           compounds_to_match_path=compounds_to_match.as_posix(),
                           output_path=self.output_folder.as_posix(),
                           tolerance=0.0015)

        ms.generate_ms_results()

        self.assertEqual(ms.tolerance, 0.0015)
        self.assertEqual(ms.compounds_to_match.shape, (269, 9))
        self.assertIsInstance(ms.matches, pd.DataFrame)
        self.assertEqual(ms.matches.shape, (4, 9))

    def test_ms_data_matcher_properties_and_setters(self):
        ms_data_path = TESTS_DATA_PATH / 'ms_data_sample' / 'ms_data.tsv'
        compounds_to_match = TESTS_DATA_PATH / 'new_compounds_sample' / 'new_compounds.tsv'
        ms = MSDataMatcher(ms_data_path=ms_data_path.as_posix(),
                           compounds_to_match_path=compounds_to_match.as_posix(),
                           output_path=self.new_output_folder.as_posix(),
                           tolerance=0.0015)

        output_path = ms.output_path
        self.assertEqual(output_path, Path(self.new_output_folder))

        ms.generate_ms_results()

        _ = ms.ms_data_path
        with self.assertRaises(FileNotFoundError):
            ms.ms_data_path = 'not_existing_path.tsv'

        ms.ms_data_path = TESTS_DATA_PATH / 'ms_data_sample' / 'ms_data_subsample.tsv'

        _ = ms.compounds_to_match
        with self.assertRaises(FileNotFoundError):
            ms.compounds_to_match = 'not_existing_path.tsv'

        ms.compounds_to_match = TESTS_DATA_PATH / 'new_compounds_sample' / 'new_compounds_subsample.tsv'

        tl = ms.tolerance
        ms.tolerance = 0.0015 + tl

        _ = ms.matches
        with self.assertRaises(AttributeError):
            ms.matches = pd.DataFrame()

