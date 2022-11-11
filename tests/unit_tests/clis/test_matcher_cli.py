import os
import shutil
from unittest import TestCase

from tests import TESTS_DATA_PATH


class MatchMSDataCLITestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        self.output_path = self.output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data.tsv')
        self.compounds_to_match_path = os.path.join(TESTS_DATA_PATH, 'results_sample/new_compounds.tsv')
        self.tolerance = 0.02

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestMatcherCLI(MatchMSDataCLITestCase, TestCase):

    def test_matcher_cli(self):
        exit_status = os.system('matcher_cli --help')
        self.assertEqual(exit_status, 0)

    def test_matcher_cli_missing_all_args(self):
        # missing argument 'MS_DATA'
        exit_status = os.system('matcher_cli')
        self.assertEqual(exit_status, 512)

    def test_matcher_cli_invalid_ms_data_path(self):
        # missing argument 'OUTPUT_PATH'
        exit_status = os.system('matcher_cli dummy_arg_1 dummy_arg_2 dummy_arg_3')
        self.assertEqual(exit_status, 256)

    def test_matcher_cli_missing_compounds_arg(self):
        # dummy argumets (FileNotFoundError)
        exit_status = os.system(f'matcher_cli {self.ms_data_path}')
        self.assertEqual(exit_status, 512)

    def test_matcher_cli_missing_output_path_arg(self):
        # dummy argumets (FileNotFoundError)
        exit_status = os.system(f'matcher_cli {self.ms_data_path} {self.compounds_to_match_path}')
        self.assertEqual(exit_status, 512)

    def test_matcher_cli_working(self):
        exit_status = os.system(f"matcher_cli {self.ms_data_path} {self.compounds_to_match_path} {self.output_path}")
        self.assertEqual(exit_status, 0)

    def test_matcher_invalid_tolerance(self):
        exit_status = os.system(f"matcher_cli {self.ms_data_path} {self.compounds_to_match_path} {self.output_path} "
                                f"--tolerance=invalid_tolerance")
        self.assertEqual(exit_status, 512)
