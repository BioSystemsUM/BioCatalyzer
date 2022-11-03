import os
from unittest import TestCase

from biocatalyzer.io_utils import Loaders

from tests import TESTS_DATA_PATH


class LoadersTestCase(TestCase):

    def test_load_compounds(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        cmps = Loaders.load_compounds(path=compounds_path)
        self.assertEqual(cmps.shape, (4, 2))
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        self.assertRaises(ValueError, Loaders.load_compounds, reaction_rules_path)

        invalid_path = 'asdasdas.tsv'
        self.assertRaises(FileNotFoundError, Loaders.load_compounds, invalid_path)

    def test_load_reaction_rules(self):
        compounds_path = os.path.join(TESTS_DATA_PATH, 'compounds_sample/compounds.tsv')
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        rules = Loaders.load_reaction_rules(path=reaction_rules_path)
        self.assertEqual(rules.shape, (51, 7))

        rules = Loaders.load_reaction_rules(path=reaction_rules_path, orgs=['eco'])
        self.assertEqual(rules.shape, (7, 7))

        self.assertRaises(ValueError, Loaders.load_reaction_rules, compounds_path)

        invalid_path = 'asdasdas.tsv'
        self.assertRaises(FileNotFoundError, Loaders.load_reaction_rules, invalid_path)

    def test_load_organisms(self):
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        organisms_path = os.path.join(TESTS_DATA_PATH, 'organisms_sample/organisms_to_use.tsv')
        orgs = Loaders.load_organisms(path=organisms_path)
        self.assertEqual(len(orgs), 2)

        self.assertRaises(ValueError, Loaders.load_organisms, reaction_rules_path)

        invalid_path = 'asdasdas.tsv'
        self.assertRaises(FileNotFoundError, Loaders.load_organisms, invalid_path)

    def test_load_byproducts_to_remove(self):
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        molecules_to_remove_path = os.path.join(TESTS_DATA_PATH, 'byproducts_to_remove_sample/byproducts.tsv')
        byproducts = Loaders.load_byproducts_to_remove(path=molecules_to_remove_path)
        self.assertEqual(len(byproducts), 9)

        self.assertRaises(ValueError, Loaders.load_byproducts_to_remove, reaction_rules_path)

        invalid_path = 'asdasdas.tsv'
        self.assertRaises(FileNotFoundError, Loaders.load_byproducts_to_remove, invalid_path)

    def test_load_patterns_to_remove(self):
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        patterns_to_remove_path = os.path.join(TESTS_DATA_PATH, 'patterns_to_remove_sample/patterns.tsv')
        patterns = Loaders.load_patterns_to_remove(path=patterns_to_remove_path)
        self.assertEqual(len(patterns), 8)

        self.assertRaises(ValueError, Loaders.load_patterns_to_remove, reaction_rules_path)

        invalid_path = 'asdasdas.tsv'
        self.assertRaises(FileNotFoundError, Loaders.load_patterns_to_remove, invalid_path)

    def test_verify_file(self):
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        self.assertTrue(Loaders._verify_file(reaction_rules_path))
        self.assertFalse(Loaders._verify_file('random_path.tsv'))

    def test_load_ms_data(self):
        ms_data_path = os.path.join(TESTS_DATA_PATH, 'ms_data_sample/ms_data.tsv')
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        self.assertRaises(FileNotFoundError, Loaders.load_ms_data, 'random_path.tsv')
        self.assertRaises(ValueError, Loaders.load_ms_data, reaction_rules_path)
        self.assertEqual(Loaders.load_ms_data(ms_data_path).shape, (438, 20))

    def test_load_new_compounds(self):
        new_compounds_path = os.path.join(TESTS_DATA_PATH, 'new_compounds_sample/new_compounds.tsv')
        reaction_rules_path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        self.assertRaises(FileNotFoundError, Loaders.load_new_compounds, 'random_path.tsv')
        self.assertRaises(ValueError, Loaders.load_new_compounds, reaction_rules_path)
        self.assertEqual(Loaders.load_new_compounds(new_compounds_path).shape, (269, 7))


