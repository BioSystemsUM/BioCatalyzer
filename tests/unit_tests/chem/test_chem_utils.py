from unittest import TestCase

from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles, MolToInchi
from rdkit.Chem.rdChemReactions import ChemicalReaction

from biocatalyzer.chem import ChemUtils


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

    def test_uncharge_mol(self):
        smiles = ['[NH3+]CC[O-]', 'CC(=O)O[IH2+2](O)OC(C)=O', '[NH3+]CC([O-])C[O-]', 'NCCO']
        for s in smiles:
            new_s = ChemUtils.uncharge_smiles(s)
            for at in MolFromSmiles(new_s).GetAtoms():
                self.assertEqual(at.GetFormalCharge(), 0)

        invalid_smiles = '[NH3+]CC[O-]]'
        self.assertIsNone(ChemUtils.uncharge_smiles(invalid_smiles))

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
