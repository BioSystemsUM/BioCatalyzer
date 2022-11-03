import os
import shutil
from unittest import TestCase

from tests import TESTS_DATA_PATH


class BioCatalyzerCLITestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        self.new_output_folder = 'new_output_path/'
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        self.output_path = self.output_folder
        self.neutralize = False
        self.reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        self.patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns.tsv')
        self.molecules_to_remove_path = os.path.join(TESTS_DATA_PATH, 'byproducts_to_remove_sample/byproducts.tsv')
        self.min_atom_count = 4
        self.n_jobs = -1
        self.ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data.tsv')
        self.match_ms_data = True
        self.tolerance = 0.02

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestBioCatalyzerCLI(BioCatalyzerCLITestCase, TestCase):

    def test_biocatalyzer_cli(self):
        exit_status = os.system('biocatalyzer_cli --help')
        self.assertEqual(exit_status, 0)

    def test_biocatalyzer_cli_missing_args(self):
        # missing argument 'COMPOUNDS'
        exit_status = os.system('biocatalyzer_cli')
        self.assertEqual(exit_status, 512)

        # missing argument 'OUTPUT_PATH'
        exit_status = os.system('biocatalyzer_cli dummy_arg_1')
        self.assertEqual(exit_status, 512)

    def test_biocatalyzer_cli_working(self):
        exit_status = os.system(f"biocatalyzer_cli {self.compounds_path} {self.output_folder}")
        self.assertEqual(exit_status, 0)

    def test_biocatalyzer_full(self):
        exit_status = os.system(f"biocatalyzer_cli {self.compounds_path} {self.output_folder} "
                                f"--neutralize={self.neutralize} --reaction_rules={self.reaction_rules_path} "
                                f"--patterns_to_remove={self.patterns_to_remove_path} "
                                f"--molecules_to_remove={self.molecules_to_remove_path} "
                                f"--min_atom_count={self.min_atom_count} --match_ms_data={self.match_ms_data} "
                                f"--ms_data_path={self.ms_data_path} --tolerance={self.tolerance} "
                                f"--n_jobs={self.n_jobs}")
        self.assertEqual(exit_status, 0)

    def test_biocatalyzer_full2(self):
        exit_status = os.system(f"biocatalyzer_cli {self.compounds_path} {self.output_folder} "
                                f"--neutralize={self.neutralize} --reaction_rules={self.reaction_rules_path} "
                                f"--patterns_to_remove={None} "
                                f"--molecules_to_remove={None} "
                                f"--min_atom_count={self.min_atom_count} --match_ms_data={self.match_ms_data} "
                                f"--ms_data_path={self.ms_data_path} --tolerance={self.tolerance} "
                                f"--n_jobs={self.n_jobs}")
        self.assertEqual(exit_status, 0)
