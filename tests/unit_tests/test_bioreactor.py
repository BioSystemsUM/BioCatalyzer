import os
import shutil
from unittest import TestCase

import pandas as pd

from biocatalyzer.bioreactor import BioReactor

from tests import TESTS_DATA_PATH


class BioReactorTestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        self.new_output_folder = 'new_output_path/'
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestBioReactor(BioReactorTestCase, TestCase):

    def test_bioreactor(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        organisms_path = os.path.join(TESTS_DATA_PATH, 'organisms_sample/organisms_to_use.tsv')
        patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns.tsv')
        molecules_to_remove_path = os.path.join(TESTS_DATA_PATH, 'byproducts_to_remove_sample/byproducts.tsv')
        br = BioReactor(compounds_path=compounds_path,
                        organisms_path=organisms_path,
                        patterns_to_remove_path=patterns_to_remove_path,
                        molecules_to_remove_path=molecules_to_remove_path,
                        output_path=self.output_folder,
                        n_jobs=12)
        br.react()

        self.assertEqual(br.reaction_rules.shape, (1368, 7))
        self.assertEqual(br.compounds.shape, (4, 2))
        self.assertIsInstance(br.new_compounds, pd.DataFrame)
        self.assertEqual(br.new_compounds.shape[1], 7)

    def test_bioreactor_all_orgs(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns.tsv')
        molecules_to_remove_path = os.path.join(TESTS_DATA_PATH, 'byproducts_to_remove_sample/byproducts.tsv')
        br_no_orgs_filter = BioReactor(compounds_path=compounds_path,
                                       patterns_to_remove_path=patterns_to_remove_path,
                                       molecules_to_remove_path=molecules_to_remove_path,
                                       output_path=self.output_folder,
                                       neutralize_compounds=True,
                                       n_jobs=12)
        br_no_orgs_filter.react()

        self.assertEqual(br_no_orgs_filter.reaction_rules.shape, (3332, 7))
        self.assertEqual(br_no_orgs_filter.compounds.shape, (4, 2))
        self.assertIsInstance(br_no_orgs_filter.new_compounds, pd.DataFrame)
        self.assertEqual(br_no_orgs_filter.new_compounds.shape[1], 7)

    def test_bioreactor_all_orgs_keep_all(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        patterns_to_remove_path = None
        molecules_to_remove_path = None
        br_no_orgs_filter = BioReactor(compounds_path=compounds_path,
                                       patterns_to_remove_path=patterns_to_remove_path,
                                       molecules_to_remove_path=molecules_to_remove_path,
                                       output_path=self.output_folder,
                                       n_jobs=-1)
        br_no_orgs_filter.react()

        self.assertEqual(br_no_orgs_filter.reaction_rules.shape, (3332, 7))
        self.assertEqual(br_no_orgs_filter.compounds.shape, (4, 2))
        self.assertIsInstance(br_no_orgs_filter.new_compounds, pd.DataFrame)
        self.assertEqual(br_no_orgs_filter.new_compounds.shape[1], 7)

    def test_bioreactor_properties_and_setters(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        organisms_path = os.path.join(TESTS_DATA_PATH, 'organisms_sample/organisms_to_use.tsv')
        br = BioReactor(compounds_path=compounds_path,
                        organisms_path=organisms_path,
                        output_path=self.output_folder,
                        n_jobs=12)

        with self.assertRaises(ValueError):
            _ = br.new_compounds
        with self.assertRaises(AttributeError):
            br.new_compounds = 'random_thing'

        output_path = br.output_path
        self.assertEqual(output_path, self.output_folder)

        br.output_path = self.new_output_folder
        shutil.rmtree(self.new_output_folder)

        with self.assertRaises(FileExistsError):
            br.output_path = os.path.join(TESTS_DATA_PATH, 'results_sample/')

        br.react()

        with self.assertRaises(FileNotFoundError):
            br.compounds = 'not_existing_path.tsv'

        br.compounds = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds_subsample.tsv')

        br.compounds = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C;C(C1C(C(C(C(O1)O)O)O)O)O'

        with self.assertRaises(FileNotFoundError):
            br.reaction_rules = 'not_existing_path.tsv'

        br.reaction_rules = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules_subsample.tsv')

        br.output_path = 'new_output_path'

        _ = br.compounds_path
        with self.assertRaises(FileNotFoundError):
            br.compounds_path = 'not_existing_path.tsv'

        br.compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds_subsample.tsv')

        _ = br.neutralize
        br.neutralize = True

        _ = br.organisms_path
        with self.assertRaises(FileNotFoundError):
            br.organisms_path = 'not_existing_path.tsv'

        br.organisms_path = os.path.join(TESTS_DATA_PATH, 'organisms_sample/organisms_subsample.tsv')

        br.organisms_path = 'hsa;eco'

        _ = br.molecules_to_remove_path
        with self.assertRaises(FileNotFoundError):
            br.molecules_to_remove_path = 'not_existing_path.tsv'

        br.molecules_to_remove_path = os.path.join(TESTS_DATA_PATH,
                                                   'byproducts_to_remove_sample/byproducts_subsample.tsv')

        _ = br.patterns_to_remove_path
        with self.assertRaises(FileNotFoundError):
            br.patterns_to_remove_path = 'not_existing_path.tsv'

        br.patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns_subsample.tsv')

        mac = br.min_atom_count
        br.min_atom_count = mac + 1

        _ = br.n_jobs
        br.n_jobs = -1
        br.n_jobs = 6

