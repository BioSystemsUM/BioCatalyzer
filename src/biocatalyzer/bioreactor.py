import logging
import multiprocessing
import os
import time
import uuid
from typing import Union

import numpy as np
import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles

from biocatalyzer.chem import ChemUtils
from biocatalyzer.io_utils import Loaders

DATA_FILES = os.path.dirname(__file__)


class BioReactor:
    """
    Main class of the BioCatalyzer package.
    Performs the transformation of the reactants into products using reaction rules.
    """

    def __init__(self,
                 compounds_path: str,
                 output_path: str,
                 neutralize_compounds: bool = False,
                 reaction_rules_path: str = 'default',
                 organisms_path: str = None,
                 molecules_to_remove_path: Union[str, None] = 'default',
                 patterns_to_remove_path: Union[str, None] = 'default',
                 min_atom_count: int = 5,
                 n_jobs: int = 1):
        """
        Initialize the BioReactor class.

        Parameters
        ----------
        compounds_path: str
            The path to the file containing the compounds to use as reactants.
        output_path: str
            The path directory to save the results to.
        neutralize_compounds: bool
            Whether to neutralize input compounds and generated compounds.
        reaction_rules_path: str
            The path to the file containing the reaction rules.
        organisms_path: str
            The path to the file containing the organisms to filter the reaction rules by.
        molecules_to_remove_path: str
            The path to the file containing the molecules to remove from the products.
        patterns_to_remove_path: str
            The path to the file containing the patterns to remove from the products.
        min_atom_count: int
            The minimum number of heavy atoms a product must have.
        n_jobs: int
            The number of jobs to run in parallel.
        """
        # silence RDKit logger
        RDLogger.DisableLog('rdApp.*')
        self._compounds_path = compounds_path
        self._output_path = output_path
        self._neutralize = neutralize_compounds
        self._organisms_path = organisms_path
        self._reaction_rules_path = reaction_rules_path
        self._molecules_to_remove_path = molecules_to_remove_path
        self._patterns_to_remove_path = patterns_to_remove_path
        self._set_up_files()
        self._orgs = Loaders.load_organisms(self._organisms_path)
        self._reaction_rules = Loaders.load_reaction_rules(self._reaction_rules_path, orgs=self._orgs)
        self._set_output_path(self._output_path)
        self._compounds = Loaders.load_compounds(self._compounds_path, self._neutralize)
        self._molecules_to_remove = Loaders.load_byproducts_to_remove(self._molecules_to_remove_path)
        self._patterns_to_remove = Loaders.load_patterns_to_remove(self._patterns_to_remove_path)
        self._min_atom_count = min_atom_count
        if n_jobs == -1:
            self._n_jobs = multiprocessing.cpu_count()
        else:
            self._n_jobs = n_jobs
        self._new_compounds = None

    @property
    def compounds(self):
        """
        Get the compounds used as reactants.

        Returns
        -------
        pd.DataFrame
            The compounds used as reactants.
        """
        return self._compounds

    @compounds.setter
    def compounds(self, compounds_path: str):
        """
        Set the compounds to use as reactants.

        Parameters
        ----------
        compounds_path: str
            The path to the file containing the compounds to use as reactants.
        """
        if compounds_path != self._compounds_path:
            self._compounds_path = compounds_path
            self._compounds = Loaders.load_compounds(self._compounds_path, self._neutralize)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def reaction_rules(self):
        """
        Get the reaction rules used.

        Returns
        -------
        pd.DataFrame
            The reaction rules used.
        """
        return self._reaction_rules

    @reaction_rules.setter
    def reaction_rules(self, reaction_rules_path: str):
        """
        Set the reaction rules to use.

        Parameters
        ----------
        reaction_rules_path: str
            The path to the file containing the reaction rules to use.
        """
        if reaction_rules_path != self._reaction_rules_path:
            self._reaction_rules = Loaders.load_reaction_rules(reaction_rules_path, orgs=self._orgs)
            self._reaction_rules_path = reaction_rules_path
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def new_compounds(self):
        """
        Get the new compounds generated.

        Returns
        -------
        pd.DataFrame
            The new compounds generated.
        """
        if isinstance(self._new_compounds, pd.DataFrame):
            return self._new_compounds
        else:
            raise ValueError('No compounds generated yet. Run the BioReactor react method first.')

    @new_compounds.setter
    def new_compounds(self, new_compounds: pd.DataFrame):
        """
        Set the new compounds generated.

        Parameters
        ----------
        new_compounds: pd.DataFrame
            The new compounds generated.
        """
        raise AttributeError('New compounds cannot be set manually! You need to run the react method!')

    @property
    def output_path(self):
        """
        Get the output path.

        Returns
        -------
        str
            The output path.
        """
        return self._output_path

    @output_path.setter
    def output_path(self, output_path: str):
        """
        Set the output path.

        Parameters
        ----------
        output_path: str
            The output path.
        """
        self._set_output_path(output_path)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def compounds_path(self):
        """
        Get the path to the compounds file.

        Returns
        -------
        str
            The path to the compounds file.
        """
        return self._compounds_path

    @compounds_path.setter
    def compounds_path(self, compounds_path: str):
        """
        Set the path to the compounds file.

        Parameters
        ----------
        compounds_path: str
            The path to the compounds file.
        """
        if compounds_path != self._compounds_path:
            self._compounds_path = compounds_path
            logging.info('Loading compounds again with the new path information...')
            self._compounds = Loaders.load_compounds(self._compounds_path, self._neutralize)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def neutralize(self):
        """
        Get whether to neutralize compounds.

        Returns
        -------
        bool
            Whether to neutralize compounds.
        """
        return self._neutralize

    @neutralize.setter
    def neutralize(self, neutralize: bool):
        """
        Set whether to neutralize compounds.

        Parameters
        ----------
        neutralize: bool
            Whether to neutralize compounds.
        """
        if neutralize != self._neutralize:
            self._neutralize = neutralize
            logging.info('Loading compounds again with the new neutralize information...')
            self._compounds = Loaders.load_compounds(self._compounds_path, self._neutralize)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def organisms_path(self):
        """
        Get the path to the organisms file.

        Returns
        -------
        str
            The path to the organisms file.
        """
        return self._organisms_path

    @organisms_path.setter
    def organisms_path(self, organisms_path: str):
        """
        Set the path to the organisms file.

        Parameters
        ----------
        organisms_path: str
            The path to the organisms file.
        """
        if organisms_path != self._organisms_path:
            self._organisms_path = organisms_path
            logging.info('Loading organisms again with the new path information...')
            self._orgs = Loaders.load_organisms(self._organisms_path)
            logging.info('Loading reaction rules again with the new organisms information...')
            self._reaction_rules = Loaders.load_reaction_rules(self._reaction_rules_path, orgs=self._orgs)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def molecules_to_remove_path(self):
        """
        Get the path to the molecules to remove file.

        Returns
        -------
        str
            The path to the molecules to remove file.
        """
        return self._molecules_to_remove_path

    @molecules_to_remove_path.setter
    def molecules_to_remove_path(self, molecules_to_remove_path: str):
        """
        Set the path to the molecules to remove file.

        Parameters
        ----------
        molecules_to_remove_path: str
            The path to the molecules to remove file.
        """
        if molecules_to_remove_path != self._molecules_to_remove_path:
            self._molecules_to_remove_path = molecules_to_remove_path
            logging.info('Loading molecules to remove again with the new path information...')
            self._molecules_to_remove = Loaders.load_byproducts_to_remove(self._molecules_to_remove_path)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def patterns_to_remove_path(self):
        """
        Get the path to the patterns to remove file.

        Returns
        -------
        str
            The path to the patterns to remove file.
        """
        return self._patterns_to_remove_path

    @patterns_to_remove_path.setter
    def patterns_to_remove_path(self, patterns_to_remove_path: str):
        """
        Set the path to the patterns to remove file.

        Parameters
        ----------
        patterns_to_remove_path: str
            The path to the patterns to remove file.
        """
        if patterns_to_remove_path != self._patterns_to_remove_path:
            self._patterns_to_remove_path = patterns_to_remove_path
            logging.info('Loading patterns to remove again with the new path information...')
            self._patterns_to_remove = Loaders.load_patterns_to_remove(self._patterns_to_remove_path)
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def min_atom_count(self):
        """
        Get the minimum atom count.

        Returns
        -------
        int
            The minimum atom count.
        """
        return self._min_atom_count

    @min_atom_count.setter
    def min_atom_count(self, min_atom_count: int):
        """
        Set the minimum atom count.

        Parameters
        ----------
        min_atom_count: int
            The minimum atom count.
        """
        if min_atom_count != self._min_atom_count:
            self._min_atom_count = min_atom_count
        if self._new_compounds is not None:
            logging.warning('Results should be generated again for the new information provided!')

    @property
    def n_jobs(self):
        """
        Get the number of jobs.

        Returns
        -------
        int
            The number of jobs.
        """
        return self._n_jobs

    @n_jobs.setter
    def n_jobs(self, n_jobs: int):
        """
        Set the number of jobs.

        Parameters
        ----------
        n_jobs: int
            The number of jobs.
        """
        if n_jobs != self._n_jobs:
            if n_jobs == -1:
                self._n_jobs = multiprocessing.cpu_count()
            else:
                self._n_jobs = n_jobs

    def _set_up_files(self):
        if self._reaction_rules_path == 'default':
            self._reaction_rules_path = os.path.join(
                DATA_FILES, 'data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv')
        if self._molecules_to_remove_path == 'default':
            self._molecules_to_remove_path = os.path.join(DATA_FILES, 'data/byproducts_to_remove/byproducts.tsv')
        if self._patterns_to_remove_path == 'default':
            self._patterns_to_remove_path = os.path.join(DATA_FILES, 'data/patterns_to_remove/patterns.tsv')

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
            if os.path.exists(output_path + '/results.tsv') or os.path.exists(output_path + '/new_compounds.tsv'):
                raise FileExistsError(f"Results in {output_path} already exists. Define a different output path so "
                                      f"that previous results are not overwritten.")

    def _match_patterns(self, smiles: str):
        """
        Check if mol matches patterns to remove.

        Parameters
        ----------
        smiles: str
            The smiles to check.

        Returns
        -------
        bool
            True if mol matches patterns to remove, False otherwise.
        """
        if len(self._patterns_to_remove) == 0:
            return False
        mol = MolFromSmiles(smiles)
        if not mol:
            return True
        for bp in self._patterns_to_remove:
            if mol.HasSubstructMatch(bp):
                return True
        return False

    def _min_atom_count_filter(self, smiles: str):
        """
        Check if mol has at least the minimum number of atoms.

        Parameters
        ----------
        smiles: str
            The smiles to check.

        Returns
        -------
        bool
            True if mol has at least the minimum number of atoms, False otherwise.
        """
        mol = MolFromSmiles(smiles)
        if not mol:
            return False
        if mol.GetNumHeavyAtoms() >= self._min_atom_count:
            return True
        else:
            return False

    def _match_byproducts(self, smiles: str):
        """
        Check if mol matches byproducts to remove.

        Parameters
        ----------
        smiles: str
            The smiles to check.

        Returns
        -------
        bool
            True if mol matches byproducts to remove, False otherwise.
        """
        if len(self._molecules_to_remove) == 0:
            return False
        if smiles in self._molecules_to_remove:
            return True
        else:
            return False

    def _get_ec_numbers(self, reaction_rule_id: str):
        """
        Get the EC numbers associated with a reaction rule.

        Parameters
        ----------
        reaction_rule_id: str
            The reaction rule id.

        Returns
        -------
        str
            The EC numbers associated with the reaction rule.
        """
        return self._reaction_rules[self._reaction_rules.InternalID == reaction_rule_id].EC_Numbers.values[0]

    @staticmethod
    def process_results(results: pd.DataFrame):
        """
        Process the results of the reactor.
        Group results by unique SMILES and merges the other columns.

        Parameters
        ----------
        results: pd.DataFrame
            The results dataframe to process.

        Returns
        -------
        pd.DataFrame
            The processed results.
        """
        results.EC_Numbers = results.EC_Numbers.fillna('')
        results = results.groupby(['OriginalCompoundID', 'NewCompoundSmiles']).agg({'OriginalCompoundSmiles': 'first',
                                                                                    'OriginalReactionRuleID': ';'.join,
                                                                                    'NewCompoundID': 'first',
                                                                                    'NewReactionSmiles': ';'.join,
                                                                                    'EC_Numbers': ';'.join}
                                                                                   ).reset_index()
        # reorder columns
        results = results[['OriginalCompoundID', 'OriginalCompoundSmiles', 'OriginalReactionRuleID', 'NewCompoundID',
                           'NewCompoundSmiles', 'NewReactionSmiles', 'EC_Numbers']]

        def merge_fields(value):
            if len(value.split(';')) == 1:
                return value
            values = []
            for v in value.split(';'):
                if v not in values:
                    values.append(v)
            return ';'.join(values)
        results['OriginalReactionRuleID'] = results['OriginalReactionRuleID'].apply(lambda x: merge_fields(x))
        results['NewReactionSmiles'] = results['NewReactionSmiles'].apply(lambda x: merge_fields(x))

        def merge_ec_numbers(x):
            if x == '':
                return np.NaN
            x = list(set(x.split(';')))
            x = [i for i in x if i != '']
            return ';'.join(x)

        results['EC_Numbers'] = results['EC_Numbers'].apply(lambda x: merge_ec_numbers(x))
        return results

    def _react_single(self, smiles: str, smarts: str):
        """
        React a single compound with a single reaction rule.
        Writes the results to the output files.

        Parameters
        ----------
        smiles: str
            The smiles of the reactant.
        smarts: str
            The SMARTS string of the reaction.
        """
        new_compounds = pd.DataFrame(columns=['OriginalCompoundID', 'OriginalCompoundSmiles', 'OriginalReactionRuleID',
                                              'NewCompoundID', 'NewCompoundSmiles', 'NewReactionSmiles', 'EC_Numbers'])
        reactants = self._reaction_rules[self._reaction_rules.SMARTS == smarts].Reactants.values[0]
        reactants = reactants.replace("Any", smiles).split(';')
        results = ChemUtils.react(reactants, smarts)
        if len(results) > 0:
            smiles_id = self._compounds[self._compounds.smiles == smiles].compound_id.values[0]
            smarts_id = self._reaction_rules[self._reaction_rules.SMARTS == smarts].InternalID.values[0]
            for i, result in enumerate(results):
                products = result.split('>')[-1].split('.')
                # keep only the most similar compound to the input compound
                most_similar_product = ChemUtils.most_similar_compound(smiles, products)
                if most_similar_product not in new_compounds.NewCompoundSmiles.values:
                    if not self._match_byproducts(most_similar_product) \
                            and not self._match_patterns(most_similar_product) \
                            and self._min_atom_count_filter(most_similar_product):
                        if self._neutralize:
                            most_similar_product = ChemUtils.uncharge_smiles(most_similar_product)
                        ecs = self._get_ec_numbers(smarts_id)
                        new_compounds.loc[len(new_compounds)] = [smiles_id, smiles, smarts_id,
                                                                 f"{smiles_id}_{uuid.uuid4()}", most_similar_product,
                                                                 result, ecs]
        return new_compounds

    def react(self):
        """
        Transform reactants into products using the reaction rules.
        """
        t0 = time.time()
        results_ = []
        for compound in self._compounds.smiles:
            with multiprocessing.Pool(self._n_jobs) as pool:
                results_.extend(pool.starmap(self._react_single, zip([compound] * self._reaction_rules.shape[0],
                                                                     self._reaction_rules.SMARTS)))

        not_empty = [not df.empty for df in results_]
        if not any(not_empty):
            logging.info('No new compounds could be generated using this reaction rules.')
            t1 = time.time()
            logging.info(f"Time elapsed: {t1 - t0} seconds")
            return False
        results = pd.concat(results_)
        results = self.process_results(results)

        results.to_csv(self._output_path + '/new_compounds.tsv', sep='\t', index=False)
        logging.info(f"New compounds saved to {self._output_path}new_compounds.tsv")
        logging.info(f"{results.shape[0]} unique new compounds generated!")
        self._new_compounds = results
        t1 = time.time()
        logging.info(f"Time elapsed: {t1 - t0} seconds")
        return True


if __name__ == '__main__':
    output_path_ = 'results/results_example/'
    br = BioReactor(compounds_path='data/compounds/drugs_paper_subset.csv',
                    output_path=output_path_,
                    organisms_path='data/organisms/organisms_to_use.tsv',
                    patterns_to_remove_path='data/patterns_to_remove/patterns.tsv',
                    molecules_to_remove_path='data/byproducts_to_remove/byproducts.tsv',
                    min_atom_count=5,
                    n_jobs=12)
    logging.basicConfig(filename=f'{output_path_}logging_bioreactor.log', level=logging.DEBUG)
    br.react()
