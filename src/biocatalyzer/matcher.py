import logging
import multiprocessing
import os
import time
from pathlib import Path
from typing import Union

import pandas as pd
from pandarallel import pandarallel

from biocatalyzer.chem import ChemUtils
from biocatalyzer.io_utils import Loaders
from biocatalyzer._utils import match_value

DATA_FILES = Path(__file__).resolve().parent


class MSDataMatcher:
    """
    Main class of the MS data matcher.
    Performs matching of the new predicted compounds to the MS data.
    """

    def __init__(self,
                 ms_data_path: str,
                 output_path: str,
                 compounds_to_match_path: str,
                 tolerance: float = 0.02,
                 n_jobs: int = 1):
        """
        Initialize the MSDataMatcher class.

        Parameters
        ----------
        ms_data_path: str
            Path to the MS data.
        output_path: str
            Path to the output directory.
        compounds_to_match_path: str
            Path to the new predicted compounds to match.
        tolerance: float
            The tolerance for the mass matching.
        n_jobs: int
            The number of jobs to run in parallel.
        """
        self._set_up_data_files(compounds_to_match_path)
        if isinstance(self._new_compounds, pd.DataFrame):
            if self._new_compounds.shape[0] == 0:
                raise ValueError('The new compounds file is empty!')
        self._ms_data_path = ms_data_path
        self._ms_data = Loaders.load_ms_data(self._ms_data_path)
        self._set_output_path(output_path)
        self._tolerance = tolerance
        if n_jobs == -1:
            self._n_jobs = multiprocessing.cpu_count()
        else:
            self._n_jobs = n_jobs
        self._calculate_masses()
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
        self._set_output_path(path)
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
            self._ms_data = Loaders.load_ms_data(self._ms_data_path)
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

    def _set_up_data_files(self, new_compounds: str):
        """
        Set up the reaction rules and new compounds data files.

        Parameters
        ----------
        new_compounds: str
            The path to the new compounds to match.
        """
        self._set_up_reaction_rules()
        self._set_up_new_compounds(new_compounds)

    def _set_up_reaction_rules(self):
        """
        Loads the reaction rules data file.
        """
        self._reaction_rules_path = DATA_FILES / 'data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv'
        self._reaction_rules = Loaders.load_reaction_rules(self._reaction_rules_path.as_posix())

    def _set_up_new_compounds(self, path: str):
        """
        Loads the new compounds data file.

        Parameters
        ----------
        path: str
            Path to the new compounds' data.
        """
        self._new_compounds = Loaders.load_new_compounds(path)

    def _set_output_path(self, output_path: str):
        """
        Make the output directory if it does not exist.

        Parameters
        ----------
        output_path: str
            The path to the output directory.
        """
        output_path = Path(output_path)
        if not output_path.exists():
            output_path.mkdir(parents=True)
        else:
            if (output_path / 'matches.tsv').exists():
                raise FileExistsError(
                    f"File {output_path / 'matches.tsv'} already exists. Define a different output path so that "
                    f"previous results are not overwritten."
                )
        self._output_path = output_path

    def _calculate_masses(self):
        """
        Calculate the masses of the new compounds.
        """
        self._new_compounds['NewCompoundExactMass'] = \
            [ChemUtils.calc_exact_mass(m) for m in self._new_compounds.NewCompoundSmiles.values]
        self._new_compounds['NewCompoundExactMass'] = self._new_compounds['NewCompoundSmiles'].apply(
            lambda x: ChemUtils.calc_exact_mass(x))

    def _match_to_parent(self, value: float, parent: str):
        """
        Match the mass value to the parent products.

        Parameters
        ----------
        value: float
            The value to match.
        parent: str
            The parent compound.

        Returns
        -------
        list:
            The list of matched indexes.
        """
        ms_data = self._ms_data[self._ms_data['ParentCompound'] == parent]
        matches = match_value(value, ms_data['Mass'].values, self._tolerance)
        return ms_data.iloc[matches].index.values

    def _match_masses(self):
        """
        Match the masses.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the matches.
        """
        ms_df = self._new_compounds
        pandarallel.initialize(nb_workers=self._n_jobs)
        ms_df['Index'] = self._new_compounds.parallel_apply(
            lambda x: self._match_to_parent(x['NewCompoundExactMass'],
                                            '_'.join(x['NewCompoundID'].split('_')[:-1])), axis=1)

        ms_df = self._new_compounds[self._new_compounds['Index'].map(lambda d: len(d)) > 0]
        ms_df = ms_df.explode('Index')
        ms_df.drop(columns=['OriginalReactionRuleID', 'NewReactionSmiles'], inplace=True)
        ms_df['ParentCompoundExactMass'] = [ChemUtils.calc_exact_mass(m) for m in ms_df.OriginalCompoundSmiles.values]
        ms_df['MassDiff'] = ms_df['ParentCompoundExactMass'] - ms_df['NewCompoundExactMass']
        ms_df = ms_df[['Index', 'OriginalCompoundID', 'OriginalCompoundSmiles', "ParentCompoundExactMass",
                       'NewCompoundID', 'NewCompoundSmiles', 'NewCompoundExactMass', 'MassDiff', 'EC_Numbers']]
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
        self._matches = self._match_masses()
        path = self._output_path / 'matches.tsv'
        self._matches.to_csv(path, sep='\t', index=False)
        logging.info(f"Matches saved to {path.as_posix()}")
        logging.info(f"{self._matches.shape[0]} matches found!")
        t1 = time.time()
        logging.info(f"Time elapsed: {t1 - t0} seconds")


if __name__ == '__main__':
    output_path_ = 'results/results_example/'
    ms = MSDataMatcher(ms_data_path='data/ms_data_example/ms_data_paper.tsv',
                       compounds_to_match_path='results/results_example/new_compounds.tsv',
                       output_path=output_path_,
                       tolerance=0.0015,
                       n_jobs=-1)
    logging.basicConfig(filename=f'{output_path_}logging_matcher.log', level=logging.DEBUG)
    ms.generate_ms_results()
