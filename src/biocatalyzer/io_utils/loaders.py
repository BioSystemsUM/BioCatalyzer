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
    def load_reaction_rules(path='data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates.tsv', orgs='ALL'):
        """
        Load the reaction rules to use.

        Parameters
        ----------
        path: str
            Path to the reaction rules.
        orgs: Union[list, str]
            List of organisms to use. If 'ALL', all organisms will be used.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the reaction rules to use.
        """
        rules = pd.read_csv(path, header=0, sep='\t')

        if 'InternalID' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "InternalID".')
        if 'Reactants' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "Reactants".')
        if 'SMARTS' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "SMARTS".')
        if 'EC_Numbers' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "EC_Numbers".')
        if 'Organisms' not in rules.columns:
            raise ValueError('The reaction rules file must contain a column named "Organisms".')

        def match_org(value, orgs_list):
            if isinstance(value, str):
                if any(s in value.split(';') for s in orgs_list):
                    return True
            return False

        if not isinstance(orgs, str):
            # TODO: check if adding expontaneous reactions actually makes sense
            orgs.append('expontaneous')
            rules['has_org'] = rules.apply(lambda x: match_org(x['Organisms'], orgs), axis=1)
            rules = rules[rules['has_org']]
            rules.drop('has_org', axis=1, inplace=True)
        return rules

    @staticmethod
    def load_organisms(path):
        """
        Load the organisms to use.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the organisms to use.
        """
        orgs = pd.read_csv(path, header=0, sep='\t')
        if 'org_id' not in orgs.columns:
            raise ValueError('The organisms file must contain a column named "org_id".')
        return orgs.org_id.values

    @staticmethod
    def load_coreactants(path='data/coreactants/all_coreactants.tsv'):
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
        if 'smiles' not in coreactants.columns:
            raise ValueError('The compounds file must contain a column named "smiles".')
        if 'compound_id' not in coreactants.columns:
            raise ValueError('The compounds file must contain a column named "compound_id".')
        return coreactants

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
