from typing import Union, List

from rdkit import Chem
from rdkit.Chem import MolFromSmiles, Mol, MolToSmiles, RemoveHs, AllChem, Descriptors
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, ChemicalReaction


class ChemUtils:
    """
    A class containing utility functions using RDKit.
    """

    @staticmethod
    def mol_to_isomerical_smiles(mol: Mol):
        """
        Converts a molecule to its canonical SMILES.

        Parameters
        ----------
        mol: Mol
            The molecule to convert.

        Returns
        -------
        str
            The SMILES string.
        """
        try:
            return MolToSmiles(RemoveHs(mol), isomericSmiles=True)
        except TypeError:
            return None

    @staticmethod
    def validate_smiles(smiles: List[str]):
        """
        Validates a list of SMILES.

        Parameters
        ----------
        smiles: List[str]
            The SMILES to validate.

        Returns
        -------
        List[str]
            The valid SMILES.
        """
        if all(MolFromSmiles(v) is None for v in smiles):
            return False
        return True

    @staticmethod
    def _smarts_to_reaction(reaction_smarts: str):
        """
        Converts a SMARTS string to a reaction.

        Parameters
        ----------
        reaction_smarts: str
            The SMARTS string.

        Returns
        -------
        ChemicalReaction
            The reaction.
        """
        try:
            return ReactionFromSmarts(reaction_smarts)
        except ValueError:
            return None

    @staticmethod
    def _remove_hs(mol: Mol):
        """
        Removes hydrogen atoms from a molecule.

        Parameters
        ----------
        mol: Mol
            The molecule to remove hydrogen atoms from.

        Returns
        -------
        Mol
            The molecule wit implicit hydrogen atoms.
        """
        try:
            return RemoveHs(mol)
        except Chem.rdchem.KekulizeException:
            return mol
        except Chem.rdchem.AtomValenceException:
            return mol
        except Chem.rdchem.AtomKekulizeException:
            return mol

    @staticmethod
    def react(smiles: Union[str, List[str]], smarts: str):
        """
        Reacts a molecule with a reaction.

        Parameters
        ----------
        smiles: Union[str, List[str]]
            The smiles of the reactant(s)' molecule(s).
        smarts: str
            The SMARTS string of the reaction.

        Returns
        -------
        list of str
            The list of products.
        """
        if isinstance(smiles, str):
            mol = (MolFromSmiles(smiles),)
        else:
            mol = [MolFromSmiles(s) for s in smiles]
        reaction = ChemUtils._smarts_to_reaction(smarts)
        if None in mol or reaction is None:
            return []
        try:
            return ChemUtils._create_reaction_instances(reaction, mol)
        except ValueError as e:
            print(e)
            return None

    @staticmethod
    def _create_reaction_instances(rxn: ChemicalReaction, reactants: List[Mol]):
        """
        Creates reaction smiles from a reaction and a list of reactants.

        Parameters
        ----------
        rxn: ChemicalReaction
            The reaction.
        reactants: List[Mol]
            The reactants.

        Returns
        -------
        list of str
            The list of reaction smiles.
        """
        res = []
        ps = rxn.RunReactants(reactants)
        for pset in ps:
            tres = ChemicalReaction()
            for p in pset:
                tres.AddProductTemplate(ChemUtils._remove_hs(p))
            for reactant in reactants:
                tres.AddReactantTemplate(ChemUtils._remove_hs(reactant))
            res.append(tres)
        return list(set([AllChem.ReactionToSmiles(entry, canonical=True) for entry in res]))

    @staticmethod
    def uncharge_smiles(smiles: str):
        """
        Neutralizes a molecule.

        Parameters
        ----------
        smiles: str
            The molecule smiles to uncharge.
        Returns
        -------
        str
            The uncharged molecule smiles.
        """
        mol = MolFromSmiles(smiles)
        if mol:
            uncharger = Uncharger()
            return MolToSmiles(uncharger.uncharge(mol))
        return None

    @staticmethod
    def calc_exact_mass(smiles: str):
        """
        Calculates the exact mass of a molecule.

        Parameters
        ----------
        smiles: str
            The molecule smiles to calculate the exact mass.
        Returns
        -------
        float
            The exact mass.
        """
        mol = MolFromSmiles(smiles)
        if mol:
            return Descriptors.ExactMolWt(mol)
        return None

    @staticmethod
    def match_masses(smiles, masses, mass_tolerance):
        mass = ChemUtils.calc_exact_mass(smiles)
        if masses is None:
            return True, mass
        if mass:
            any_found = any(mass - mass_tolerance <= m <= mass + mass_tolerance for m in masses)
            return any_found, round(mass, 4)
        return False, mass
