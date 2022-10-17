import os
from typing import Union

import pandas as pd

from biocatalyzer.chem import ChemUtils
from biocatalyzer.io_utils import Loaders
from biocatalyzer._utils import match_value

DATA_FILES = os.path.dirname(__file__)


class MSDataMatcher:

    def __init__(self,
                 ms_data_path: str,
                 output_path: str,
                 compounds_to_match: Union[pd.DataFrame, str],
                 compound_id_field: str,
                 compound_smiles_field: str,
                 mode: str = 'mass',
                 ms_field: str = 'Mass',
                 tolerance: float = 0.02):
        self._ms_data = Loaders.load_ms_data(os.path.join(DATA_FILES, ms_data_path), ms_field)
        self._output_path = output_path
        self._compound_id_field = compound_id_field
        self._compound_smiles_field = compound_smiles_field
        self._set_output_path(self._output_path)
        self._mode = mode
        self._ms_field = ms_field
        self._tolerance = tolerance
        self._set_up_data_files(compounds_to_match)
        self._prepare_mode()

    def _prepare_mode(self):
        """
        Processes the new compounds data according to the mode.
        """
        if self._mode == 'mass':
            self._new_compounds['NewCompoundExactMass'] = \
                [ChemUtils.calc_exact_mass(m) for m in self._new_compounds.NewCompoundSmiles.values]
        elif self._mode == 'mass_diff':
            self._new_compounds['NewCompoundExactMass'] = \
                [ChemUtils.calc_exact_mass(m) for m in self._new_compounds.NewCompoundSmiles.values]
            self._new_compounds['NewCompoundExactMassDiff'] = \
                [ChemUtils.calc_exact_mass(row['OriginalCompoundSmiles']) -
                 ChemUtils.calc_exact_mass(row['NewCompoundSmiles']) for i, row in self._new_compounds.iterrows()]
        else:
            raise ValueError('The mode must be either "mass" or "mass_dif".')

    def _set_up_data_files(self, new_compounds: Union[pd.DataFrame, str]):
        """
        Set up the reaction rules and new compounds data files.

        Parameters
        ----------
        new_compounds: Union[pd.DataFrame, str]
            The new compounds to match.
        """
        self._set_up_reaction_rules()
        if not isinstance(new_compounds, pd.DataFrame):
            self._set_up_new_compounds(new_compounds)
        else:
            self._new_compounds = new_compounds

    def _set_up_reaction_rules(self):
        """
        Loads the reaction rules data file.
        """
        self._reaction_rules_path = \
            os.path.join(DATA_FILES, 'data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates.tsv')
        self._reaction_rules = Loaders.load_reaction_rules(self._reaction_rules_path)

    def _set_up_new_compounds(self, path: str):
        """
        Loads the new compounds data file.

        Parameters
        ----------
        path: str
            Path to the new compounds' data.
        """
        self._new_compounds = Loaders.load_new_compounds(path)

    @staticmethod
    def _set_output_path(output_path: str):
        """
        Make the output directory if it does not exist.

        Parameters
        ----------
        output_path: str
            The path to the output directory.
        """
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        else:
            if os.path.exists(output_path + '/matches.tsv'):
                raise FileExistsError(f"File {output_path} already exists. Define a different output path so that "
                                      f"previous results are not overwritten.")

    def _match_masses(self):
        """
        Match the masses.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the matches.
        """
        ms_df = pd.DataFrame(columns=[self._compound_id_field, self._compound_smiles_field,
                                      f"{self._compound_id_field}_ExactMass", self._ms_field, 'NewCompoundID',
                                      'NewCompoundSmiles', 'NewCompoundExactMass', 'EC_Numbers'])
        for i, row in self._new_compounds.iterrows():
            mv, mi = match_value(row['NewCompoundExactMass'], self._ms_data[self._ms_field].values, self._tolerance)
            if mv and self._ms_data.loc[mi, self._compound_id_field] == row['NewCompoundID'].split('_')[0]:
                ms_df.loc[len(ms_df)] = [self._ms_data.loc[mi, self._compound_id_field],
                                         self._ms_data.loc[mi, self._compound_smiles_field],
                                         ChemUtils.calc_exact_mass(self._ms_data.loc[mi, self._compound_smiles_field]),
                                         self._ms_data.loc[mi, self._ms_field],
                                         row['NewCompoundID'],
                                         row['NewCompoundSmiles'],
                                         row['NewCompoundExactMass'],
                                         row['EC_Numbers']]
        return ms_df

    def _match_mass_diff(self):
        """
        Match mass differences.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the matches.
        """
        ms_df = pd.DataFrame(columns=[self._compound_id_field, self._compound_smiles_field,
                                      f"{self._compound_id_field}_ExactMass", self._ms_field, 'NewCompoundID',
                                      'NewCompoundSmiles', 'NewCompoundExactMass', 'NewCompoundExactMassDiff',
                                      'EC_Numbers'])
        for i, row in self._new_compounds.iterrows():
            mv, mi = match_value(row['NewCompoundExactMassDiff'], self._ms_data[self._ms_field].values, self._tolerance)
            if mv and self._ms_data.loc[mi, self._compound_id_field] == row['NewCompoundID'].split('_')[0]:
                ms_df.loc[len(ms_df)] = [self._ms_data.loc[mi, self._compound_id_field],
                                         self._ms_data.loc[mi, self._compound_smiles_field],
                                         ChemUtils.calc_exact_mass(self._ms_data.loc[mi, self._compound_smiles_field]),
                                         self._ms_data.loc[mi, self._ms_field],
                                         row['NewCompoundID'],
                                         row['NewCompoundSmiles'],
                                         row['NewCompoundExactMass'],
                                         row['NewCompoundExactMassDiff'],
                                         row['EC_Numbers']]
        return ms_df

    def generate_ms_results(self):
        """
        Get the dataset with the matches.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the matches.
        """
        if self._mode == 'mass':
            matches = self._match_masses()
        elif self._mode == 'mass_diff':
            matches = self._match_mass_diff()
        else:
            raise ValueError('The mode must be either "mass" or "mass_dif".')
        matches.to_csv(self._output_path + '/matches.tsv', sep='\t', index=False)


if __name__ == '__main__':
    ms = MSDataMatcher(ms_data_path='data/ms_data_example/ms_data_paper.tsv',
                       compounds_to_match='results/results_example/new_compounds.tsv',
                       output_path='results/results_example',
                       compound_id_field='ParentDrug',
                       compound_smiles_field='SMILES',
                       mode='mass',
                       ms_field='MZ',
                       tolerance=0.0015)

    ms.generate_ms_results()
