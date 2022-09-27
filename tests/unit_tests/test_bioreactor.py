import os
import shutil
from unittest import TestCase

from mock_bioreactor import MockBioReactor


class BioReactorTestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        os.mkdir(self.output_folder)

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestBioReactor(BioReactorTestCase, TestCase):

    def test_bioreactor(self):
        br = MockBioReactor(compounds_path='data/compounds_sample/compounds.tsv',
                            reaction_rules_path='data/reaction_rules_sample/reactionrules.tsv',
                            coreactants_path='data/coreactants_sample/coreactants.tsv',
                            patterns_to_remove_path='data/patterns_to_remove_sample/patterns.tsv',
                            molecules_to_remove_path='data/byproducts_to_remove_sample/byproducts.tsv',
                            output_path=self.output_folder,
                            n_jobs=12)
        br.react()
