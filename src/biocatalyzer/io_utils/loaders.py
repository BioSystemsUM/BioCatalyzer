import os

import pandas as pd
from rdkit.Chem import MolFromSmarts, MolFromSmiles

from biocatalyzer.chem import ChemUtils


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
            Path to the compounds or string with ;-separated SMILES.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the compounds to use.
        """

        if Loaders._verify_file(path):
            compounds = pd.read_csv(path, header=0, sep='\t')
            if 'smiles' not in compounds.columns:
                raise ValueError('The compounds file must contain a column named "smiles".')
            if 'compound_id' not in compounds.columns:
                raise ValueError('The compounds file must contain a column named "compound_id".')
            return compounds[['compound_id', 'smiles']]
        elif ChemUtils.validate_smiles(path.split(';')):
            df = pd.DataFrame()
            df['smiles'] = path.split(';')
            df['compound_id'] = [f'input_compound_{i}' for i in range(len(df))]
            return df
        else:
            raise FileNotFoundError(f"File {path} not found.")

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

        Parameters
        ----------
        path: str
            Path to the organisms or ;-separated list of organisms identifiers.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the organisms to use.
        """
        if not path:
            return 'ALL'
        if Loaders._verify_file(path):
            orgs = pd.read_csv(path, header=0, sep='\t')
            if 'org_id' not in orgs.columns:
                raise ValueError('The organisms file must contain a column named "org_id".')
            return list(orgs.org_id.values)
        else:
            return path.split(';')

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
        if not path:
            return []
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
        if not path:
            return []
        patterns = pd.read_csv(path, header=0, sep='\t')
        if 'smarts' not in patterns.columns:
            raise ValueError('The coreactants file must contain a column named "smarts".')
        return [MolFromSmarts(sp) for sp in patterns.smarts.values]

    @staticmethod
    def load_masses_to_match(masses):
        """
        Load masses to match.

        Parameters
        ----------
        masses: str
            Path to the masses to match.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the masses to match.
        """
        if not masses:
            return None
        elif Loaders._verify_file(masses):
            masses = pd.read_csv(masses, header=0, sep='\t')
            if 'mass' not in masses.columns:
                raise ValueError('The masses file must contain a column named "mass".')
            return list(masses['mass'].values)
        else:
            return [float(m) for m in masses.split(';')]

    @staticmethod
    def _verify_file(path: str):
        """
        Verify that the provided paths to the files exist.

        Parameters
        ----------
        path: str
            The path to the files to verify.
        """
        if not os.path.exists(path):
            return False
        return True
