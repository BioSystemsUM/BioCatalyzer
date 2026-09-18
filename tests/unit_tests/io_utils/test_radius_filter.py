import os
import tempfile
from unittest import TestCase

import pandas as pd

from biocatalyzer.io_utils import Loaders

from tests import TESTS_DATA_PATH


class RadiusFilterTestCase(TestCase):
    """
    Tests for the radius filter of Loaders.load_reaction_rules.

    RetroRules v3.0.0 models each template at one or more radii between 0 and 10, listed
    in a `Radii` field. The filter must handle the full range, and a rule belongs to a
    radius when that radius appears in its `Radii` list: one SMARTS pattern can be
    generated at several radii, so radius-specific subsets are not disjoint.

    The fixture is built in the test rather than added to tests/data, so the expected
    membership of every radius is defined alongside the assertions.
    """

    RADII = list(range(11))

    def setUp(self):
        """
        Build a rule set where the membership of each radius is known by construction.

        Rule i is modelled at every radius that is a multiple of (i + 1), so:
            rule 0 -> all radii 0..10
            rule 1 -> 0, 2, 4, 6, 8, 10
            rule 2 -> 0, 3, 6, 9
            ...
        Radius 0 therefore matches all eleven rules, and each radius has a membership
        that can be computed independently of the implementation.
        """
        self._tmp = tempfile.TemporaryDirectory()
        rows = []
        for i in range(11):
            radii = [r for r in self.RADII if r % (i + 1) == 0]
            rows.append({
                'InternalID': f'Rule_{i}',
                'Reactants': 'Any',
                'SMARTS': f'[C:{i + 1}]>>[C:{i + 1}]',
                'EC_Numbers': '1.1.1.1',
                'Organisms': 'eco;bsu',
                'Radii': ','.join(str(r) for r in radii),
            })
        self.expected = {
            r: {f'Rule_{i}' for i in range(11) if r % (i + 1) == 0}
            for r in self.RADII
        }
        self.rules_path = os.path.join(self._tmp.name, 'rules_with_radii.tsv')
        pd.DataFrame(rows).to_csv(self.rules_path, sep='\t', index=False)

    def tearDown(self):
        self._tmp.cleanup()

    def test_every_radius_from_0_to_10(self):
        """Each of the eleven radii selects exactly the rules modelled at it."""
        for radius in self.RADII:
            rules = Loaders.load_reaction_rules(path=self.rules_path, radius=radius)
            self.assertEqual(set(rules['InternalID']), self.expected[radius],
                             msg=f'wrong selection for radius {radius}')

    def test_radius_accepts_int_and_string(self):
        """An integer and its string form select the same rules."""
        for radius in self.RADII:
            as_int = Loaders.load_reaction_rules(path=self.rules_path, radius=radius)
            as_str = Loaders.load_reaction_rules(path=self.rules_path, radius=str(radius))
            self.assertEqual(set(as_int['InternalID']), set(as_str['InternalID']))

    def test_radius_list(self):
        """A ;-separated list selects the union of the individual radii."""
        rules = Loaders.load_reaction_rules(path=self.rules_path, radius='4;6;8')
        expected = self.expected[4] | self.expected[6] | self.expected[8]
        self.assertEqual(set(rules['InternalID']), expected)

    def test_radius_range(self):
        """A range is inclusive at both ends."""
        rules = Loaders.load_reaction_rules(path=self.rules_path, radius='4:6')
        expected = self.expected[4] | self.expected[5] | self.expected[6]
        self.assertEqual(set(rules['InternalID']), expected)

        full = Loaders.load_reaction_rules(path=self.rules_path, radius='0:10')
        self.assertEqual(len(full), 11)

    def test_subsets_are_not_disjoint(self):
        """
        A template modelled at several radii belongs to all of them.

        This is the property that makes membership, rather than equality with a single
        modelled radius, the correct test: the sizes of the radius subsets sum to more
        than the number of rules.
        """
        total = sum(len(Loaders.load_reaction_rules(path=self.rules_path, radius=r))
                    for r in self.RADII)
        self.assertGreater(total, 11)

    def test_default_is_radius_6(self):
        """
        The default radius is 6, not 'ALL'.

        Radius 6 was adopted as the operating point: on the case study it removes 99.4%
        of the products generated at radius 0 while keeping 30 of the 33 input compounds
        productive. Calling the loader without a radius must therefore select the same
        rules as asking for radius 6 explicitly, and not the whole rule set.
        """
        implicit = Loaders.load_reaction_rules(path=self.rules_path)
        explicit = Loaders.load_reaction_rules(path=self.rules_path, radius=6)
        self.assertEqual(set(implicit['InternalID']), set(explicit['InternalID']))
        self.assertEqual(set(implicit['InternalID']), self.expected[6])

    def test_all_disables_the_filter(self):
        """
        'ALL' is how a caller opts out of the default.

        This needed no test while the default was 'ALL'; it does now, because it is the
        only way to recover the unfiltered rule set.
        """
        rules = Loaders.load_reaction_rules(path=self.rules_path, radius='ALL')
        self.assertEqual(len(rules), 11)

    def test_radius_composes_with_organisms(self):
        """The radius and organism filters apply together."""
        rules = Loaders.load_reaction_rules(path=self.rules_path, orgs=['eco'], radius=5)
        self.assertEqual(set(rules['InternalID']), self.expected[5])

        rules = Loaders.load_reaction_rules(path=self.rules_path, orgs=['hsa'], radius=5)
        self.assertEqual(len(rules), 0)

    def test_rule_set_without_radii_column_is_untouched(self):
        """
        Rule sets that declare no radii are returned unfiltered, with a warning.

        This is what keeps the filter backwards compatible with the rule files
        distributed with the tool, which carry no `Radii` column.
        """
        path = os.path.join(TESTS_DATA_PATH, 'reaction_rules_sample/reactionrules.tsv')
        unfiltered = Loaders.load_reaction_rules(path=path)
        for radius in self.RADII:
            with self.assertLogs(level='WARNING'):
                rules = Loaders.load_reaction_rules(path=path, radius=radius)
            self.assertEqual(len(rules), len(unfiltered))

    def test_parse_radius(self):
        """The radius specification is parsed from every accepted form."""
        self.assertEqual(Loaders._parse_radius(6), [6])
        self.assertEqual(Loaders._parse_radius('6'), [6])
        self.assertEqual(Loaders._parse_radius('4;6;8'), [4, 6, 8])
        self.assertEqual(Loaders._parse_radius('4:8'), [4, 5, 6, 7, 8])
        self.assertEqual(Loaders._parse_radius('0:10'), list(range(11)))
        self.assertEqual(Loaders._parse_radius([2, 4]), [2, 4])
