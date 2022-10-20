import os
import shutil
from unittest import TestCase

from tests import TESTS_DATA_PATH


class BioReactorCLITestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        os.mkdir(self.output_folder)

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestBioReactorCLI(BioReactorCLITestCase, TestCase):

    def test_bioreactor_cli(self):
        exit_status = os.system('bioreactor_cli --help')
        self.assertEqual(exit_status, 0)

    def test_bioreactor_cli_missing_args(self):
        # missing argument 'COMPOUNDS'
        exit_status = os.system('bioreactor_cli')
        self.assertEqual(exit_status, 512)

        # missing argument 'OUTPUT_PATH'
        exit_status = os.system('bioreactor_cli dummy_arg_1')
        self.assertEqual(exit_status, 512)

    def test_bioreactor_cli_dummy_args(self):
        # dummy argumets (FileNotFoundError)
        exit_status = os.system('bioreactor_cli dummy_arg_1 dummy_arg_2')
        self.assertEqual(exit_status, 256)

    def test_bioreactor_cli_working(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        exit_status = os.system(f"bioreactor_cli {compounds_path} {self.output_folder}")
        self.assertEqual(exit_status, 0)

    def test_bioreactor_invalid_args(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        output_path = self.output_folder
        neutralize = False
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        organisms_path = os.path.join(TESTS_DATA_PATH, 'organisms_sample/organisms_to_use.tsv')
        patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns.tsv')
        molecules_to_remove_path = os.path.join(TESTS_DATA_PATH, 'byproducts_to_remove_sample/byproducts.tsv')
        min_atom_count = 4
        n_jobs = -1

        cli = f"bioreactor_cli {compounds_path} {output_path} --neutralize={neutralize} " \
              f"--reaction_rules={reaction_rules_path} --organisms={organisms_path} " \
              f"--patterns_to_remove={patterns_to_remove_path} --molecules_to_remove={molecules_to_remove_path} " \
              f"--min_atom_count={min_atom_count} --n_jobs={n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, 0)
        shutil.rmtree(self.output_folder)

        # invalid neutralize value
        neutralize = 10
        cli = f"bioreactor_cli {compounds_path} {output_path} --neutralize={neutralize} " \
              f"--reaction_rules={reaction_rules_path} --organisms={organisms_path} " \
              f"--patterns_to_remove={patterns_to_remove_path} --molecules_to_remove={molecules_to_remove_path} " \
              f"--min_atom_count={min_atom_count} --n_jobs={n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, 512)

