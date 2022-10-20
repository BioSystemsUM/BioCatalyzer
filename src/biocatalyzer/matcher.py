import logging
import os
import time
from typing import Union

import pandas as pd

from biocatalyzer.chem import ChemUtils
from biocatalyzer.io_utils import Loaders
from biocatalyzer._utils import match_value

DATA_FILES = os.path.dirname(__file__)


class MSDataMatcher:
    """
    Main class of the MS data matcher.
    Performs matching of the new predicted compounds to the MS data.
    """

    def __init__(self,
                 ms_data_path: str,
                 output_path: str,
                 compounds_to_match: Union[pd.DataFrame, str],
                 mode: str = 'mass',
                 tolerance: float = 0.02):
        """
        Initialize the MSDataMatcher class.

        Parameters
        ----------
        ms_data_path: str
            Path to the MS data.
        output_path: str
            Path to the output directory.
        compounds_to_match: Union[pd.DataFrame, str]
            The new predicted compounds to match.
        mode: str
            The mode of the matcher. Either 'mass' or 'mass_diff'.
        tolerance: float
            The tolerance for the mass matching.
        """
        self._ms_data_path = ms_data_path
        self._ms_data = Loaders.load_ms_data(self._ms_data_path, mode)
        self._output_path = output_path
        self._set_output_path(self._output_path)
        self._mode = mode
        self._tolerance = tolerance
        self._set_up_data_files(compounds_to_match)
        self._prepare_mode()
        self._matches = None

    @property
    def output_path(self):
        """
        Returns the output path.

        Returns
        -------
        str:
            The output path.
        """
        return self._output_path

    @output_path.setter
    def output_path(self, path: str):
        """
        Set the output path.

        Parameters
        ----------
        path: str
            The output path.
        """
        self._output_path = path
        self._set_output_path(self._output_path)
        if self._matches is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def ms_data_path(self):
        """
        Returns the path to the MS data.

        Returns
        -------
        str:
            The path to the MS data.
        """
        return self._ms_data_path

    @ms_data_path.setter
    def ms_data_path(self, path: str):
        """
        Set the path to the MS data.

        Parameters
        ----------
        path: str
            The path to the MS data.
        """
        if path != self._ms_data_path:
            self._ms_data_path = path
            logging.info('Loading MS data with the new path information...')
            self._ms_data = Loaders.load_ms_data(self._ms_data_path, self._mode)
        if self._matches is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def compounds_to_match(self):
        """
        Returns the compounds to match.

        Returns
        -------
        pd.DataFrame:
            The compounds to match.
        """
        return self._new_compounds

    @compounds_to_match.setter
    def compounds_to_match(self, new_compounds: Union[pd.DataFrame, str]):
        """
        Set the compounds to match.

        Parameters
        ----------
        new_compounds: Union[pd.DataFrame, str]
            The new compounds to match.
        """
        self._set_up_data_files(new_compounds)
        logging.info('Loading the new compounds to match with the new path information...')
        self._prepare_mode()
        if self._matches is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def mode(self):
        """
        Returns the mode of the matcher.

        Returns
        -------
        str:
            The mode of the matcher.
        """
        return self._mode

    @mode.setter
    def mode(self, mode: str):
        """
        Set the mode of the matcher.

        Parameters
        ----------
        mode: str
            The mode of the matcher.
        """
        if mode != self._mode:
            self._mode = mode
            self._ms_data = Loaders.load_ms_data(self._ms_data_path, mode)
            self._prepare_mode()
        if self._matches is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def tolerance(self):
        """
        Returns the tolerance.

        Returns
        -------
        float:
            The tolerance.
        """
        return self._tolerance

    @tolerance.setter
    def tolerance(self, tolerance: float):
        """
        Set the tolerance.

        Parameters
        ----------
        tolerance: float
            The tolerance.
        """
        if tolerance != self._tolerance:
            self._tolerance = tolerance
        if self._matches is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def matches(self):
        """
        Returns the matches.

        Returns
        -------
        pd.DataFrame:
            The matches.
        """
        return self._matches

    @matches.setter
    def matches(self, matches: pd.DataFrame):
        """
        Set the matches.

        Parameters
        ----------
        matches: pd.DataFrame
            The matches.
        """
        raise AttributeError('Matches cannot be set manually! You need to run the generate_ms_results method!')

    def _prepare_mode(self):
        """
        Processes the new compounds data according to the mode.
        """
        if self._mode == 'mass':
            self._ms_field = 'Mass'
            self._new_compounds['NewCompoundExactMass'] = \
                [ChemUtils.calc_exact_mass(m) for m in self._new_compounds.NewCompoundSmiles.values]
        elif self._mode == 'mass_diff':
            self._ms_field = 'MassDiff'
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
        self._reaction_rules_path = os.path.join(
            DATA_FILES, 'data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv')
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
        ms_df = pd.DataFrame(columns=['ParentCompound', 'ParentCompoundSmiles', "ParentCompound_ExactMass",
                                      self._ms_field, 'NewCompoundID', 'NewCompoundSmiles', 'NewCompoundExactMass',
                                      'EC_Numbers'])
        for i, row in self._new_compounds.iterrows():
            mv, mi = match_value(row['NewCompoundExactMass'], self._ms_data[self._ms_field].values, self._tolerance)
            if mv and self._ms_data.loc[mi, 'ParentCompound'] == '_'.join(row['NewCompoundID'].split('_')[:-1]):
                ms_df.loc[len(ms_df)] = [self._ms_data.loc[mi, 'ParentCompound'],
                                         self._ms_data.loc[mi, 'ParentCompoundSmiles'],
                                         ChemUtils.calc_exact_mass(self._ms_data.loc[mi, 'ParentCompoundSmiles']),
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
        ms_df = pd.DataFrame(columns=['ParentCompound', 'ParentCompoundSmiles', "ParentCompound_ExactMass",
                                      self._ms_field, 'NewCompoundID', 'NewCompoundSmiles', 'NewCompoundExactMass',
                                      'NewCompoundExactMassDiff', 'EC_Numbers'])
        for i, row in self._new_compounds.iterrows():
            mv, mi = match_value(row['NewCompoundExactMassDiff'], self._ms_data[self._ms_field].values, self._tolerance)
            if mv and self._ms_data.loc[mi, 'ParentCompound'] == row['NewCompoundID'].split('_')[0]:
                ms_df.loc[len(ms_df)] = [self._ms_data.loc[mi, 'ParentCompound'],
                                         self._ms_data.loc[mi, 'ParentCompoundSmiles'],
                                         ChemUtils.calc_exact_mass(self._ms_data.loc[mi, 'ParentCompoundSmiles']),
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
        t0 = time.time()
        if self._mode == 'mass':
            self._matches = self._match_masses()
        elif self._mode == 'mass_diff':
            self._matches = self._match_mass_diff()
        else:
            raise ValueError('The mode must be either "mass" or "mass_dif".')
        self._matches.to_csv(self._output_path + '/matches.tsv', sep='\t', index=False)
        logging.info(f"Matches saved to {self._output_path}/matches.tsv")
        logging.info(f"{self._matches.shape[0]} matches found!")
        t1 = time.time()
        logging.info(f"Time elapsed: {t1 - t0} seconds")


if __name__ == '__main__':
    output_path_ = 'results/results_example/'
    ms = MSDataMatcher(ms_data_path='data/ms_data_example/ms_data_paper.tsv',
                       compounds_to_match='results/results_example/new_compounds.tsv',
                       output_path=output_path_,
                       mode='mass',
                       tolerance=0.0015)
    logging.basicConfig(filename=f'{output_path_}logging_matcher.log', level=logging.DEBUG)
    ms.generate_ms_results()
