import os
import shutil
from unittest import TestCase

from biocatalyzer.bioreactor import BioReactor


class BioReactorTestCase(TestCase):

    def setUp(self):
        self.output_folder = 'results/'
        os.mkdir(self.output_folder)

    def tearDown(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)


class TestBioReactor(BioReactorTestCase, TestCase):

    def test_bioreactor(self):
        br = BioReactor(compounds_path='tests/data/compounds_sample/compounds.tsv',
                        reaction_rules_path='tests/data/reaction_rules_sample/reactionrules.tsv',
                        organisms_path='tests/data/organisms_sample/organisms_to_use.tsv',
                        coreactants_path='tests/data/coreactants_sample/coreactants.tsv',
                        patterns_to_remove_path='tests/data/patterns_to_remove_sample/patterns.tsv',
                        molecules_to_remove_path='tests/data/byproducts_to_remove_sample/byproducts.tsv',
                        output_path=self.output_folder,
                        n_jobs=12)
        br.react()

    def test_bioreactor_all_orgs(self):
        br_no_orgs_filter = BioReactor(compounds_path='tests/data/compounds_sample/compounds.tsv',
                                       reaction_rules_path='tests/data/reaction_rules_sample/reactionrules.tsv',
                                       coreactants_path='tests/data/coreactants_sample/coreactants.tsv',
                                       patterns_to_remove_path='tests/data/patterns_to_remove_sample/patterns.tsv',
                                       molecules_to_remove_path='tests/data/byproducts_to_remove_sample/byproducts.tsv',
                                       output_path=self.output_folder,
                                       n_jobs=12)
        br_no_orgs_filter.react()
