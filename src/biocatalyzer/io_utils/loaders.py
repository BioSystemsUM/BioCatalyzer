import logging
import os
from pathlib import Path
from typing import Union, List

import pandas as pd
from rdkit.Chem import MolFromSmarts, MolFromSmiles

from biocatalyzer.chem import ChemUtils


class Loaders:
    """
    Class containing a set of input utilities.
    """

    @staticmethod
    def load_compounds(path: str, neutralize: bool = False):
        """
        Load compounds to use.

        Parameters
        ----------
        path: str
            Path to the compounds or string with ;-separated SMILES.
        neutralize: bool
            Whether to neutralize compounds or not.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the compounds to use.
        """
        if Loaders._verify_file(path):
            path = Path(path)
            compounds = pd.read_csv(path, header=0, sep='\t')
            if 'smiles' not in compounds.columns:
                raise ValueError('The compounds file must contain a column named "smiles".')
            if 'compound_id' not in compounds.columns:
                raise ValueError('The compounds file must contain a column named "compound_id".')
            return compounds[['compound_id', 'smiles']]
        elif ChemUtils.validate_smiles(path.split(';')):
            df = pd.DataFrame()
            df['smiles'] = path.split(';')
            if neutralize:
                df['smiles'] = [ChemUtils.uncharge_smiles(s) for s in df.smiles.values]
            df['compound_id'] = [f'input_compound_{i}' for i in range(len(df))]
            return df
        else:
            raise FileNotFoundError(f"File {path} not found.")

    @staticmethod
    def load_reaction_rules(path: str,
                            orgs: Union[str, List[str]] = 'ALL',
                            radius: Union[str, int, List[int]] = 'ALL') -> pd.DataFrame:
        """
        Load the reaction rules to use.

        Parameters
        ----------
        path: str
            Path to the reaction rules.
        orgs: Union[list, str]
            List of organisms to use. If 'ALL', all organisms will be used.
        radius: Union[str, int, List[int]]
            Reaction rule radius (or radii) to keep. If 'ALL', no radius filter is applied.
            Accepts an integer (6), a ;-separated list ('4;6;8') or a range ('4:8').
            A rule is kept when any of these radii appears in its 'Radii' field, the
            same membership test already used for 'Organisms'. Only applied when the
            reaction rules file provides a 'Radii' column; rule sets without it
            (e.g. the bundled ones) are left untouched.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the reaction rules to use.
        """
        if not Loaders._verify_file(path):
            raise FileNotFoundError(f"File {path} not found.")
        path = Path(path)
        if path.suffix == '.bz2':
            rules = pd.read_csv(path, header=0, sep='\t', compression='bz2')
        else:
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
            orgs.append('spontaneous_reaction')
            rules['has_org'] = rules.apply(lambda x: match_org(x['Organisms'], orgs), axis=1)
            rules = rules[rules['has_org']]
            rules.drop('has_org', axis=1, inplace=True)

        if not (isinstance(radius, str) and radius == 'ALL'):
            if 'Radii' not in rules.columns:
                logging.warning("A radius filter was requested but the reaction rules file has no "
                                "'Radii' column. The filter was ignored.")
            else:
                radii = Loaders._parse_radius(radius)

                def match_radius(value, radii_list):
                    if isinstance(value, str):
                        return any(int(r) in radii_list for r in value.split(',') if r != '')
                    return False

                rules = rules[rules['Radii'].apply(lambda v: match_radius(v, radii))]
                logging.info(f'Using {len(rules)} reaction rules modelled at radius in {sorted(radii)}.')
        return rules

    @staticmethod
    def _parse_radius(radius: Union[str, int, List[int]]) -> List[int]:
        """
        Parse the radius specification into an explicit list of radii.

        Parameters
        ----------
        radius: Union[str, int, List[int]]
            An integer (6), a ;-separated list ('4;6;8') or an inclusive range ('4:8').

        Returns
        -------
        List[int]:
            The radii to keep.
        """
        if isinstance(radius, int):
            return [radius]
        if isinstance(radius, (list, tuple)):
            return [int(r) for r in radius]
        radius = str(radius).strip()
        if ':' in radius:
            start, end = radius.split(':')
            return list(range(int(start), int(end) + 1))
        return [int(r) for r in radius.split(';') if r != '']

    @staticmethod
    def load_organisms(path: str) -> Union[str, List[str]]:
        """
        Load the organisms to use.

        Parameters
        ----------
        path: str
            Path to the organisms or ;-separated list of organisms identifiers.

        Returns
        -------
        Union[str, List[str]]:
            List of organisms identifiers.
        """
        if path is None or path == 'None':
            return 'ALL'
        if Loaders._verify_file(path):
            path = Path(path)
            orgs = pd.read_csv(path, header=0, sep='\t')
            if 'org_id' not in orgs.columns:
                raise ValueError('The organisms file must contain a column named "org_id".')
            logging.info(f'Using {orgs.org_id.to_list()} as the Organisms.')
            return orgs.org_id.to_list()
        elif len(path.split('.')) > 1:
            raise FileNotFoundError(f"File {path} not found.")
        else:
            logging.info(f'Using {path.split(";")} as the Organisms.')
            return path.split(';')

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
        if path is None or path == 'None':
            return []
        path = Path(path)
        byproducts = pd.read_csv(path, header=0, sep='\t')
        if 'smiles' not in byproducts.columns:
            raise ValueError('The molecules to remove file must contain a column named "smiles".')
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
        if path is None or path == 'None':
            return []
        path = Path(path)
        patterns = pd.read_csv(path, header=0, sep='\t')
        if 'smarts' not in patterns.columns:
            raise ValueError('The patterns to remove file must contain a column named "smarts".')
        return [MolFromSmarts(sp) for sp in patterns.smarts.values]

    @staticmethod
    def _verify_file(path: str):
        """
        Verify that the provided paths to the files exist.

        Parameters
        ----------
        path: str
            The path to the files to verify.

        Returns
        -------
        bool:
            True if the file exists, False otherwise.
        """
        if not Path(path).exists():
            return False
        return True

    @staticmethod
    def load_ms_data(path: str):
        """
        Load the MS data.

        Parameters
        ----------
        path: str
            Path to the MS data.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the MS data.
        """
        if Loaders._verify_file(path):
            ms_data = pd.read_csv(path, header=0, sep='\t')
            if 'ParentCompound' not in ms_data.columns:
                raise ValueError('The MS data file must contain a column named "ParentCompound".')
            if 'ParentCompoundSmiles' not in ms_data.columns:
                raise ValueError('The MS data file must contain a column named "ParentCompoundSmiles".')
            if 'Mass' not in ms_data.columns:
                raise ValueError('The MS data file must contain a column named "Mass".')
            return ms_data
        else:
            raise FileNotFoundError(f"File {path} not found.")

    @staticmethod
    def load_new_compounds(path: str):
        """
        Load the new compounds data to match with the MS data.
        The file must be a new_compounds.tsv file resulting from running the BioCatalyzer BioReactor.

        Parameters
        ----------
        path: str
            Path to the new compounds data.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the new compounds' data.
        """
        if Loaders._verify_file(path):
            new_compounds = pd.read_csv(path, header=0, sep='\t')
            columns = ['OriginalCompoundID', 'OriginalCompoundSmiles', 'OriginalReactionRuleID', 'NewCompoundID',
                       'NewCompoundSmiles', 'NewReactionSmiles', 'EC_Numbers']
            if not all(col in new_compounds.columns for col in columns):
                raise ValueError(f'The new compounds file must be a result of BioCatalyzer module, i.e. it should '
                                 f'contain the following columns: {columns}.')
            return new_compounds
        else:
            raise FileNotFoundError(f"File {path} not found.")
