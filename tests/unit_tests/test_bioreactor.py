import os
import shutil
from unittest import TestCase

from biocatalyzer.bioreactor import BioReactor

from tests import TESTS_DATA_PATH


class BioReactorTestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        os.mkdir(self.output_folder)

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

    def test_bioreactor_all_orgs_no_remove_anything(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        patterns_to_remove_path = None
        molecules_to_remove_path = None
        br_no_orgs_filter = BioReactor(compounds_path=compounds_path,
                                       patterns_to_remove_path=patterns_to_remove_path,
                                       molecules_to_remove_path=molecules_to_remove_path,
                                       output_path=self.output_folder,
                                       n_jobs=12)
        br_no_orgs_filter.react()
