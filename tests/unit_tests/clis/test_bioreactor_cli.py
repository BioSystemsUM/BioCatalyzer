import os
import platform
import shutil
from unittest import TestCase

from tests import TESTS_DATA_PATH


class BioReactorCLITestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        self.output_path = self.output_folder
        self.neutralize = False
        self.reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        self.organisms_path = os.path.join(TESTS_DATA_PATH, 'organisms_sample/organisms_to_use.tsv')
        self.patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns.tsv')
        self.molecules_to_remove_path = os.path.join(TESTS_DATA_PATH, 'byproducts_to_remove_sample/byproducts.tsv')
        self.min_atom_count = 4
        self.n_jobs = -1

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestBioReactorCLI(BioReactorCLITestCase, TestCase):

    def test_bioreactor_cli(self):
        exit_status = os.system('bioreactor_cli --help')
        self.assertEqual(exit_status, 0)

    def test_bioreactor_cli_missing_args(self):
        expected_exit_code = 512 if platform.system() != 'Windows' else 2
        # missing argument 'COMPOUNDS'
        exit_status = os.system('bioreactor_cli')
        self.assertEqual(exit_status, expected_exit_code)

        # missing argument 'OUTPUT_PATH'
        exit_status = os.system('bioreactor_cli dummy_arg_1')
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_cli_dummy_args(self):
        expected_exit_code = 256 if platform.system() != 'Windows' else 1
        # dummy argumets (FileNotFoundError)
        exit_status = os.system('bioreactor_cli dummy_arg_1 dummy_arg_2')
        self.assertEqual(exit_status, expected_exit_code)
        shutil.rmtree('dummy_arg_2')

    def test_bioreactor_cli_working(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        exit_status = os.system(f"bioreactor_cli {compounds_path} {self.output_folder}")
        self.assertEqual(exit_status, 0)

    def test_bioreactor_valid_args(self):
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, 0)

    def test_bioreactor_compounds_string(self):
        compounds = "CC=C(=O)CCC(=O)O;COC(=O)C(C)CC;CCCCCC"
        cli = f"bioreactor_cli '{compounds}' {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"
        exit_status = os.system(cli)
        self.assertEqual(exit_status, 0)

    def test_bioreactor_invalid_neutralize(self):
        expected_exit_code = 512 if platform.system() != 'Windows' else 2
        # invalid neutralize value
        invalid_neutralize = 10
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={invalid_neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_invalid_reaction_rules_path(self):
        expected_exit_code = 256 if platform.system() != 'Windows' else 1
        # invalid reaction rules path
        invalid_reaction_rules_path = 'random_path'
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={invalid_reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_invalid_organisms_path(self):
        expected_exit_code = 256 if platform.system() != 'Windows' else 1
        # invalid organisms' path (it will use only spontaneous reactions)
        invalid_organisms_path = 'random_path'
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={invalid_organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, 0)

        # invalid organisms' path (it will recognize the string as a path, and it will raise a FileNotFoundError)
        invalid_organisms_path = 'random_string_as_a_path.tsv'
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={invalid_organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_valid_organisms_string(self):
        # valid organisms but in string format
        string_organisms = "hsa;eco"
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms='{string_organisms}' " \
              f"--patterns_to_remove={self.patterns_to_remove_path} --molecules_to_remove={self.molecules_to_remove_path} " \
              f"--min_atom_count={self.min_atom_count} --n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, 0)

    def test_bioreactor_invalid_patterns_to_remove_path(self):
        expected_exit_code = 256 if platform.system() != 'Windows' else 1
        # invalid patterns to remove path
        invalid_patterns_to_remove_path = 'random_path'
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={invalid_patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_invalid_molecules_to_remove_path(self):
        expected_exit_code = 256 if platform.system() != 'Windows' else 1
        # invalid molecules to remove path
        invalid_molecules_to_remove_path = 'random_path'
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={invalid_molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_invalid_min_atom_count(self):
        expected_exit_code = 512 if platform.system() != 'Windows' else 2
        # invalid min atom count
        invalid_min_atom_count = True
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={invalid_min_atom_count} " \
              f"--n_jobs={self.n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)

    def test_bioreactor_invalid_n_jobs(self):
        expected_exit_code = 512 if platform.system() != 'Windows' else 2
        # invalid n_jobs
        invalid_n_jobs = True
        cli = f"bioreactor_cli {self.compounds_path} {self.output_path} --neutralize={self.neutralize} " \
              f"--reaction_rules={self.reaction_rules_path} --organisms={self.organisms_path} " \
              f"--patterns_to_remove={self.patterns_to_remove_path} " \
              f"--molecules_to_remove={self.molecules_to_remove_path} --min_atom_count={self.min_atom_count} " \
              f"--n_jobs={invalid_n_jobs}"

        exit_status = os.system(cli)
        self.assertEqual(exit_status, expected_exit_code)
