from unittest import TestCase

from biocatalyzer.io_utils import Loaders


class LoadersTestCase(TestCase):

    def test_load_compounds(self):
        cmps = Loaders.load_compounds(path='data/compounds_sample/compounds.tsv')
        self.assertEqual(cmps.shape, (3, 2))

        self.assertRaises(ValueError, Loaders.load_compounds, 'data/reaction_rules_sample/reactionrules.tsv')

    def test_load_reaction_rules(self):
        rules = Loaders.load_reaction_rules(path='data/reaction_rules_sample/reactionrules.tsv')
        self.assertEqual(rules.shape, (51, 7))

        rules = Loaders.load_reaction_rules(path='data/reaction_rules_sample/reactionrules.tsv', orgs=['eco'])
        self.assertEqual(rules.shape, (30, 7))

        self.assertRaises(ValueError, Loaders.load_reaction_rules, 'data/compounds_sample/compounds.tsv')

    def test_load_organisms(self):
        orgs = Loaders.load_organisms(path='data/organisms_sample/organisms_to_use.tsv')
        self.assertEqual(len(orgs), 2)

        self.assertRaises(ValueError, Loaders.load_organisms, 'data/reaction_rules_sample/reactionrules.tsv')

    def test_load_coreactants(self):
        coreactants = Loaders.load_coreactants(path='data/coreactants_sample/coreactants.tsv')
        self.assertEqual(coreactants.shape, (16, 2))

        self.assertRaises(ValueError, Loaders.load_coreactants, 'data/reaction_rules_sample/reactionrules.tsv')

    def test_load_byproducts_to_remove(self):
        byproducts = Loaders.load_byproducts_to_remove(path='data/byproducts_to_remove_sample/byproducts.tsv')
        self.assertEqual(len(byproducts), 9)

        self.assertRaises(ValueError, Loaders.load_byproducts_to_remove, 'data/reaction_rules_sample/reactionrules.tsv')

    def test_load_patterns_to_remove(self):
        patterns = Loaders.load_patterns_to_remove(path='data/patterns_to_remove_sample/patterns.tsv')
        self.assertEqual(len(patterns), 8)

        self.assertRaises(ValueError, Loaders.load_patterns_to_remove, 'data/reaction_rules_sample/reactionrules.tsv')
