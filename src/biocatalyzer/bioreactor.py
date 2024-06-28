import itertools
import logging
import multiprocessing
import os
import time
import uuid
from typing import Union

import pandas as pd
from tqdm import tqdm

from biocatalyzer._utils import _merge_fields
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
        ChemUtils.rdkit_logs(False)
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
        self._new_compounds_path = os.path.join(self._output_path, 'new_compounds.tsv')
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
        if self._new_compounds is not None:
            return Loaders.load_compounds(self._new_compounds_path, False)
        else:
            raise ValueError('No compounds generated yet. Run the BioReactor react method first.')

    @new_compounds.setter
    def new_compounds(self, new_compounds: str):
        """
        Set the new compounds generated.

        Parameters
        ----------
        new_compounds: str
            The path to the file containing the new compounds generated.
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
                DATA_FILES, 'data/reactionrules/reaction_rules_biocatalyzer.tsv.bz2')
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
        mol = ChemUtils.smiles_to_mol(smiles)
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
        mol = ChemUtils.smiles_to_mol(smiles)
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

    def _match_conditions(self, smiles: str):
        """
        Check if mol matches conditions to remove.

        Parameters
        ----------
        smiles: str
            The smiles to check.

        Returns
        -------
        bool
            True if mol matches conditions to remove, False otherwise.
        """
        if not smiles:
            return False
        if '*' in smiles:
            return False
        if self._min_atom_count > 0:
            if not self._min_atom_count_filter(smiles):
                return False
        if len(self._molecules_to_remove) > 0:
            if self._match_byproducts(smiles):
                return False
        if len(self._patterns_to_remove) > 0:
            if self._match_patterns(smiles):
                return False
        return True

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

    def process_results(self, save: bool = True, overwrite: bool = True):
        """
        Process the results of the reactor.
        Group results by unique SMILES and merges the other columns.

        Parameters
        ----------
        save: bool
            If True, save the results to a file.
        overwrite: bool
            If True, overwrite the results file if it already exists.

        Returns
        -------
        Tuple[pd.DataFrame, str]
            The processed results and the path to the results file.
        """
        results = pd.read_csv(self._new_compounds_path, sep='\t', header=0)
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

        results['OriginalReactionRuleID'] = results['OriginalReactionRuleID'].apply(lambda x: _merge_fields(x))
        results['NewReactionSmiles'] = results['NewReactionSmiles'].apply(lambda x: _merge_fields(x))
        results['EC_Numbers'] = results['EC_Numbers'].apply(lambda x: _merge_fields(x))
        if save:
            if overwrite:
                results_file_proc = os.path.join(self._output_path, 'new_compounds.tsv')
                results.to_csv(results_file_proc, sep='\t', index=False)
            else:
                results_file_proc = os.path.join(self._output_path, 'new_compounds_processed.tsv')
                results.to_csv(results_file_proc, sep='\t', index=False)
        else:
            results_file_proc = self._new_compounds_path
        return results, results_file_proc

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
        reactants = self._reaction_rules[self._reaction_rules.SMARTS == smarts].Reactants.values[0]
        reactants = reactants.replace("Any", smiles).split(';')
        results = ChemUtils.react(reactants, smarts)
        if len(results) > 0:
            smiles_id = self._compounds[self._compounds.smiles == smiles].compound_id.values[0]
            smarts_id = self._reaction_rules[self._reaction_rules.SMARTS == smarts].InternalID.values[0]
            most_similar_products_set = set()
            for i, result in enumerate(results):
                products = result.split('>')[-1].split('.')
                # keep only the most similar compound to the input compound
                most_similar_product = ChemUtils.most_similar_compound(smiles, products)
                most_similar_product = ChemUtils.smiles_to_isomerical_smiles(most_similar_product)
                if most_similar_product not in most_similar_products_set:
                    most_similar_products_set.add(most_similar_product)
                    if self._match_conditions(most_similar_product):
                        if self._neutralize:
                            most_similar_product = ChemUtils.uncharge_smiles(most_similar_product)
                        ecs = self._get_ec_numbers(smarts_id)
                        with open(self._new_compounds_path, 'a',  newline='', encoding='utf-8') as f:
                            f.write(f"{smiles_id}\t{smiles}\t{smarts_id}\t{smiles_id}_{uuid.uuid4()}\t"
                                    f"{most_similar_product}\t{result}\t{ecs}\n")
                            if f"{smiles_id}\t{smiles}\t{smarts_id}\t{smiles_id}_{uuid.uuid4()}\t{most_similar_product}\t{result}\t{ecs}\n".split('\t') != 7:
                                e = f"{smiles_id}\t{smiles}\t{smarts_id}\t{smiles_id}_{uuid.uuid4()}\t{most_similar_product}\t{result}\t{ecs}\n"
                                raise ValueError('Wrong number of columns. Got: ', e)


    def react(self):
        """
        Transform reactants into products using the reaction rules.
        """
        t0 = time.time()
        with open(self._new_compounds_path, 'w') as f:
            f.write('OriginalCompoundID\tOriginalCompoundSmiles\tOriginalReactionRuleID\tNewCompoundID\t'
                    'NewCompoundSmiles\tNewReactionSmiles\tEC_Numbers\n')
        params = list(itertools.product(self._compounds.smiles, self._reaction_rules.SMARTS))
        with multiprocessing.Pool(self._n_jobs) as pool:
            pool.starmap(self._react_single, tqdm(params, total=len(params)))
        self._new_compounds = f"New products saved to {self._new_compounds_path}"
        t1 = time.time()
        logging.info(f"Time elapsed: {t1 - t0} seconds")


if __name__ == '__main__':
    output_path_ = 'results/results_example/'
    br = BioReactor(compounds_path='data/compounds/drugs.csv',
                    output_path=output_path_,
                    organisms_path='data/organisms/organisms_to_use.tsv',
                    patterns_to_remove_path='data/patterns_to_remove/patterns.tsv',
                    molecules_to_remove_path='data/byproducts_to_remove/byproducts.tsv',
                    min_atom_count=5,
                    n_jobs=12)
    logging.basicConfig(filename=f'{output_path_}logging_bioreactor.log', level=logging.DEBUG)
    br.react()
