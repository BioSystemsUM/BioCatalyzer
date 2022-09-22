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
    def load_reaction_rules(path):
        """
        Load the reaction rules to use.

        Parameters
        ----------
        path: str
            Path to the reaction rules.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the reaction rules to use.
        """
        rules = pd.read_csv(path, header=0, sep='\t')
        if 'rule_id' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "rule_id".')
        if 'smarts' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "smarts".')
        if 'coreactants_ids' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "coreactants_ids".')
        if 'reaction_info' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "reaction_info".')
        return rules[['rule_id', 'smarts', 'coreactants_ids', 'reaction_info']]

    @staticmethod
    def load_coreactants(path):
        """
        Load the coreactants to use.

        Parameters
        ----------
        path: str
            Path to the coreactants.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the coreactants to use.
        """
        coreactants = pd.read_csv(path, header=0, sep='\t')
        if 'compound_id' not in coreactants.columns:
            raise ValueError('The coreactants file must contain a column named "compound_id".')
        if 'smiles' not in coreactants.columns:
            raise ValueError('The coreactants file must contain a column named "smiles".')
        return coreactants[['compound_id', 'smiles']]

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
