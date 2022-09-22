from typing import Union, List

from rdkit import Chem
from rdkit.Chem import MolFromSmiles, Mol, MolToSmiles, RemoveHs, AllChem
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
        Removes hydrogens from a molecule.

        Parameters
        ----------
        mol: Mol
            The molecule to remove hydrogens from.

        Returns
        -------
        Mol
            The molecule wit implicit hydrogens.
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
        return ChemUtils._create_reaction_instances(reaction, mol)

    @staticmethod
    def _create_reaction_instances(rxn, reactants):
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
