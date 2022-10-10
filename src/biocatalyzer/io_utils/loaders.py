import pandas as pd
from rdkit.Chem import MolFromSmarts, MolFromSmiles


class Loaders:
    """
    Class containing a set of input utilities.
    """

    @staticmethod
    def load_compounds(path):
        """
        Load compounds to use.

        Parameters
        ----------
        path: str
            Path to the compounds.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the compounds to use.
        """
        compounds = pd.read_csv(path, header=0, sep='\t')
        if 'smiles' not in compounds.columns:
            raise ValueError('The compounds file must contain a column named "smiles".')
        if 'compound_id' not in compounds.columns:
            raise ValueError('The compounds file must contain a column named "compound_id".')
        return compounds[['compound_id', 'smiles']]

    @staticmethod
    def load_reaction_rules():
        """
        Load the reaction rules to use.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the reaction rules to use.
        """
        return pd.read_csv('data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates.tsv', header=0, sep='\t')

    @staticmethod
    def load_coreactants():
        """
        Load the coreactants to use.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the coreactants to use.
        """
        return pd.read_csv('data/coreactants/all_coreactants.tsv', header=0, sep='\t')

    @staticmethod
    def load_byproducts_to_remove(path):
        """
        Load byproducts to remove from products.

        Parameters
        ----------
        path: str
            Path to the byproducts to remove.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the byproducts to remove.
        """
        byproducts = pd.read_csv(path, header=0, sep='\t')
        if 'smiles' not in byproducts.columns:
            raise ValueError('The coreactants file must contain a column named "smiles".')
        return [MolFromSmiles(sp) for sp in byproducts.smiles.values]

    @staticmethod
    def load_patterns_to_remove(path):
        """
        Load patterns to remove from products.

        Parameters
        ----------
        path: str
            Path to the patterns to remove.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the patterns to remove.
        """
        patterns = pd.read_csv(path, header=0, sep='\t')
        if 'smarts' not in patterns.columns:
            raise ValueError('The coreactants file must contain a column named "smarts".')
        return [MolFromSmarts(sp) for sp in patterns.smarts.values]
