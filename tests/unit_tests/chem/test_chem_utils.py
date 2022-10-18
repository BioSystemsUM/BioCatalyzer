from unittest import TestCase

from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles, MolToInchi
from rdkit.Chem.rdChemReactions import ChemicalReaction

from biocatalyzer.chem import ChemUtils
from biocatalyzer.chem._utils import _correct_number_of_parenthesis


class TestChemUtils(TestCase):
    # mute rdkit logs
    RDLogger.DisableLog('rdApp.*')

    def test_mol_to_isomerical_smiles(self):
        smiles = ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                  'C(C1C(C(C(C(O1)O)O)O)O)O',
                  'CC(=O)OC1=CC=CC=C1C(=O)O']
        mols = [MolFromSmiles(s) for s in smiles]
        invalid_smiles = 'C(C1C(C(C(C(O1)O)O)O)O)O('
        invalid_mol = MolFromSmiles(invalid_smiles)

        def same_compound(mol1, mol2):
            return MolToInchi(mol1) == MolToInchi(mol2)

        for i, m in enumerate(mols):
            self.assertNotEqual(ChemUtils.mol_to_isomerical_smiles(m), smiles[i])
            self.assertTrue(same_compound(m, MolFromSmiles(ChemUtils.mol_to_isomerical_smiles(m))))

        self.assertIsNone(ChemUtils.mol_to_isomerical_smiles(invalid_mol))

    def test_validate_smiles(self):
        smiles = ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                  'C(C1C(C(C(C(O1)O)O)O)O)O',
                  'CC(=O)OC1=CC=CC=C1C(=O)O(']
        self.assertTrue(ChemUtils.validate_smiles(smiles))
        self.assertFalse(ChemUtils.validate_smiles([smiles[2]]))

    def test_smarts_to_reaction(self):
        smarts = ['[#6:10]-[#7H2,#16H1:9].[O:7]=[#6:3]-1-[#6;h1:4]=[#6:5]-[#6:6]=[#6:1]-[#6:2]-1=[O:8]>>[#6:10]-[*:9]-[#6:4]-1=[#6:5]-[#6:6]=[#6:1]-[#6:2](-[#8:8])=[#6:3]-1-[#8:7]',
                  '[#7H2,SH1:11]-[#6:10]-[#6:9]-[#6:6]-1=[#6:1]-[#6:2](=[O:8])-[#6:3](=[O:7])-[#6:4]=[#6;h1:5]-1>>[#8:8]-[#6:2]-1=[#6:3](-[#8:7])-[#6:4]=[#6:5]-2-[*:11]-[#6:10]-[#6:9]-[#6:6]-2=[#6:1]-1',
                  '[#6:1]-[#6H1:2]=[O:3].[#8:4]-[#8:5]>>[#6:1]-[#6:2](-[#8:5])=[O:3].[#8:4]']
        reactions = [ChemUtils._smarts_to_reaction(s) for s in smarts]
        invalid_smarts = '[#6:1]-[#6H1:2]=[O:3].[#8:4]-[#8:5]>>[#6:1]-[#6:2](-[#8:5])=[O:3].[#8:4]['
        invalid_reaction = ChemUtils._smarts_to_reaction(invalid_smarts)

        for r in reactions:
            self.assertIsInstance(r, ChemicalReaction)

        self.assertIsNone(invalid_reaction)

    def test_remove_hs(self):
        smiles = ['[H]c1oc2c(=O)c([H])c([H])c(=O)c=2oc1[H]',
                  '[H]c1nc([H])c([H])c([H])c1[H]']

        def check_if_molecule_has_hydrogens(mol):
            for atom in mol.GetAtoms():
                atomic_num = atom.GetAtomicNum()
                if atomic_num == 1:
                    return True
            return False

        for s in smiles:
            self.assertFalse(check_if_molecule_has_hydrogens(ChemUtils._remove_hs(MolFromSmiles(s))))

    def test_react(self):
        smiles = ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                  'C(C1C(C(C(C(O1)O)O)O)O)O',
                  'CC(=O)OC1=CC=CC=C1C(=O)O']
        smarts = [
            '[#6:10]-[#7H2,#16H1:9].[O:7]=[#6:3]-1-[#6;h1:4]=[#6:5]-[#6:6]=[#6:1]-[#6:2]-1=[O:8]>>[#6:10]-[*:9]-[#6:4]-1=[#6:5]-[#6:6]=[#6:1]-[#6:2](-[#8:8])=[#6:3]-1-[#8:7]',
            '[#7H2,SH1:11]-[#6:10]-[#6:9]-[#6:6]-1=[#6:1]-[#6:2](=[O:8])-[#6:3](=[O:7])-[#6:4]=[#6;h1:5]-1>>[#8:8]-[#6:2]-1=[#6:3](-[#8:7])-[#6:4]=[#6:5]-2-[*:11]-[#6:10]-[#6:9]-[#6:6]-2=[#6:1]-1',
            '[#6:1]-[#6H1:2]=[O:3].[#8:4]-[#8:5]>>[#6:1]-[#6:2](-[#8:5])=[O:3].[#8:4]']

        for s in smarts:
            self.assertIsNone(ChemUtils.react(smiles, s))

        known_reactant = 'Nc1nc(NC2CC2)c2ncn(C3C=CC(CO)C3)c2n1'
        coreactant = 'O=C1C=CC=CC1=O'
        known_rule = '[#6:10]-[#7H2,#16H1:9].[O:7]=[#6:3]-1-[#6;h1:4]=[#6:5]-[#6:6]=[#6:1]-[#6:2]-1=[O:8]>>[#6:10]-[*:9]-[#6:4]-1=[#6:5]-[#6:6]=[#6:1]-[#6:2](-[#8:8])=[#6:3]-1-[#8:7]'

        reaction_smiles = ChemUtils.react([known_reactant, coreactant], known_rule)
        self.assertEqual(known_reactant, reaction_smiles[0].split('.')[0])
        self.assertEqual(coreactant, reaction_smiles[0].split('.')[1].split('>>')[0])

    def test_create_reaction_instances(self):
        rule_smarts = '[#7:1].[#8:2].[#8:3]=[#6:4]1-[#6:5]=[#6:6]-[#6:7](=[#8:8])-[#6:9]=[#6:10]-1>>[#7+:1]-[#8:2].[#8:3]-[#6:4]1:[#6:10]:[#6:9]:[#6:7](-[#8:8]):[#6:6]:[#6:5]:1'
        rxn = ChemUtils._smarts_to_reaction(rule_smarts)
        reactants_smiles = 'C[NH+](C)CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12;O;*C1=C(*)C(=O)C(*)=C(*)C1=O'
        reactants = [MolFromSmiles(s) for s in reactants_smiles.split(';')]
        reaction_instances = ChemUtils._create_reaction_instances(rxn, reactants)
        self.assertEqual(3, len(reaction_instances))
        for instance in reaction_instances:
            self.assertEqual(3, len(instance.split('>')))
            for reac in reactants_smiles.split(';'):
                self.assertTrue(reac in instance.split('>')[0])

    def test_uncharge_smiles(self):
        smiles = ['[NH3+]CC[O-]', 'CC(=O)O[IH2+2](O)OC(C)=O', '[NH3+]CC([O-])C[O-]', 'NCCO']
        for s in smiles:
            new_s = ChemUtils.uncharge_smiles(s)
            for at in MolFromSmiles(new_s).GetAtoms():
                self.assertEqual(at.GetFormalCharge(), 0)

        invalid_smiles = '[NH3+]CC[O-]]'
        self.assertEqual(ChemUtils.uncharge_smiles(invalid_smiles), '[NH3+]CC[O-]]')

    def test_calc_exact_mass(self):
        smiles = ['[NH3+]CC[O-]', 'CC(=O)O[IH2+2](O)OC(C)=O', '[NH3+]CC([O-])C[O-]', 'NCCO']
        for s in smiles:
            self.assertIsInstance(ChemUtils.calc_exact_mass(s), float)

        invalid_smiles = '[NH3+]CC[O-]]'
        self.assertIsNone(ChemUtils.calc_exact_mass(invalid_smiles))

    def test_match_mass(self):
        smiles = ['CC(=O)O[IH2+2](O)OC(C)=O', '[NH3+]CC[O-]', 'NCCO', '[NH3+]CC([O-])C[O-]']
        masses = [ChemUtils.calc_exact_mass(s) for s in smiles]

        for s in smiles:
            self.assertTrue(ChemUtils.match_masses(s, masses, mass_tolerance=0))

        self.assertFalse(ChemUtils.match_masses(smiles[0], masses[1:], mass_tolerance=0.02)[0])

        masses = [263.9684]
        self.assertTrue(ChemUtils.match_masses(smiles[0], masses, mass_tolerance=0.02)[0])
        self.assertFalse(ChemUtils.match_masses(smiles[0], masses, mass_tolerance=0.01)[0])

    def test_calc_fingerprint_similarity(self):
        smile = 'C[NH+](C)CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12'
        smiles_to_compare = ['CO',
                             'C[NH2+]CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12',
                             'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1OP(=O)(O)O',
                             'CC=(']
        sims = [ChemUtils.calc_fingerprint_similarity(smile, stc) for stc in smiles_to_compare]
        self.assertEqual(4, len(sims))
        self.assertEqual(0.0, sims[-1])
        self.assertGreater(sims[1], sims[0])
        self.assertGreater(sims[1], sims[2])
        self.assertGreater(sims[1], sims[3])

    def test_most_similar_compound(self):
        smile = 'C[NH+](C)CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12'
        smiles_to_compare = ['CO',
                             'C[NH2+]CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12',
                             'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1OP(=O)(O)O']

        self.assertEqual(ChemUtils.most_similar_compound(smile, smiles_to_compare), smiles_to_compare[1])

        smiles_to_compare.append('CC=(')
        self.assertEqual(ChemUtils.most_similar_compound(smile, smiles_to_compare), smiles_to_compare[1])

        smiles_to_compare_2 = ["*c1c(*)c(O)c(*)c(*)c1O", "C[N+](C)(O)CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12"]

        self.assertEqual(ChemUtils.most_similar_compound(smile, smiles_to_compare_2), smiles_to_compare_2[1])

        smiles2 = 'CCCC(=O)Nc1ccc(OCC(O)C[NH2+]C(C)C)c(C(C)=O)c1'
        # mol number 2 is invalid
        smiles_to_compare_3 = ['*c1c(*)c(O)c(*)c(*)c1O', 'CCCC(=O)Nc1ccc(=OCC(O)C[NH2+]C(C)C)c(C(C)=O)c1']

        self.assertEqual(ChemUtils.most_similar_compound(smiles2, smiles_to_compare_3), smiles_to_compare_3[0])

    def test_correct_number_of_parenthesis(self):
        smiles = ['CO)',
                  'C[NH2+]CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12',
                  '(Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1OP(=O)(O)O']
        corrected_smiles = _correct_number_of_parenthesis(smiles)
        self.assertEqual(corrected_smiles[0], 'CO')
        self.assertEqual(corrected_smiles[1], 'C[NH2+]CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12')
        self.assertEqual(corrected_smiles[2],
                         'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1OP(=O)(O)O')
