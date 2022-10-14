import multiprocessing
import os
import time
from typing import List, Union

import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles

from biocatalyzer.chem import ChemUtils
from biocatalyzer.io_utils import Loaders


class BioReactor:
    """
    Main class of the BioCatalyzer package.
    Performs the transformation of the reactants into products using reaction rules.
    """

    def __init__(self,
                 compounds_path: str,
                 output_path: str,
                 neutralize_compounds: bool = False,
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
        self._molecules_to_remove_path = molecules_to_remove_path
        self._patterns_to_remove_path = patterns_to_remove_path
        self._set_up_files()
        self._orgs = Loaders.load_organisms(self._organisms_path)
        self._reaction_rules = Loaders.load_reaction_rules(self._reaction_rules_path, orgs=self._orgs)
        self._set_output_path(self._output_path)
        self._compounds = Loaders.load_compounds(self._compounds_path, neutralize_compounds)
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

    @property
    def new_compounds(self):
        """
        Get the new compounds generated.

        Returns
        -------
        pd.DataFrame
            The new compounds generated.
        """
        if self._new_compounds:
            return self._new_compounds
        else:
            raise ValueError('No compounds generated yet. Run the BioReactor react method first.')

    def _set_up_files(self):
        self.DATA_FILES = os.path.join(os.path.dirname(__file__), 'data')
        self._reaction_rules_path = \
            os.path.join(self.DATA_FILES, 'reactionrules/all_reaction_rules_forward_no_smarts_duplicates.tsv')
        if self._molecules_to_remove_path == 'default':
            self._molecules_to_remove_path = os.path.join(self.DATA_FILES, 'byproducts_to_remove/byproducts.tsv')
        if self._patterns_to_remove_path == 'default':
            self._patterns_to_remove_path = os.path.join(self.DATA_FILES, 'patterns_to_remove/patterns.tsv')

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
                raise FileExistsError(f"File {output_path} already exists. Define a different output path so that "
                                      f"previous results are not overwritten.")

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
        new_compounds = pd.DataFrame(columns=['ResultID', 'OriginalCompoundID', 'OriginalReactionRuleID',
                                              'NewCompoundID', 'NewCompoundSMILES', 'EC_Numbers'])
        reactants = self._reaction_rules[self._reaction_rules.SMARTS == smarts].Reactants.values[0]
        reactants = reactants.replace("Any", smiles).split(';')
        results = ChemUtils.react(reactants, smarts)
        if len(results) > 0:
            smiles_id = self._compounds[self._compounds.smiles == smiles].compound_id.values[0]
            smarts_id = self._reaction_rules[self._reaction_rules.SMARTS == smarts].InternalID.values[0]
            with open(self._output_path + '/results.tsv', 'a') as rs:
                for i, result in enumerate(results):
                    id_result = f"{smiles_id}_{smarts_id}_{i}"
                    new_product_generated = False
                    products = result.split('>')[-1].split('.')
                    for p in products:
                        # deal with cases where invalid number of parentheses are generated
                        if (p.count('(') + p.count(')')) % 2 != 0:
                            if p[0] == '(':
                                p = p[1:]
                            elif p[-1] == ')':
                                p = p[:-1]
                        if p not in new_compounds.NewCompoundSMILES.values:
                            if not self._match_byproducts(p) and not self._match_patterns(p) \
                                    and self._min_atom_count_filter(p):
                                if self._neutralize:
                                    p = ChemUtils.uncharge_smiles(p)
                                ecs = self._get_ec_numbers(smarts_id)
                                new_compounds.loc[len(new_compounds)] = [id_result, smiles_id, smarts_id,
                                                                         f"{id_result}_{i}", p, ecs]
                                new_product_generated = True
                    if new_product_generated:
                        rs.write(f"{id_result}\t{smiles_id}\t{smarts_id}\t{result}\n")
        return new_compounds

    def react(self):
        """
        Transform reactants into products using the reaction rules.
        """
        t0 = time.time()
        with open(self._output_path + '/results.tsv', 'w') as rs:
            rs.write('ResultID\tCompoundID\tRuleID\tNewReactionSmiles\n')
        for compound in self._compounds.smiles:
            with multiprocessing.Pool(self._n_jobs) as pool:
                results_ = pool.starmap(self._react_single, zip([compound]*self._reaction_rules.shape[0],
                                                                self._reaction_rules.SMARTS))

        results = pd.concat(results_)
        results = results.drop_duplicates(subset=['NewCompoundSMILES'])
        results.to_csv(self._output_path + '/new_compounds.tsv', sep='\t', index=False)
        print(f"{results.shape[0]} unique new compounds generated!")
        self._new_compounds = results
        t1 = time.time()
        print(f"Time elapsed: {t1 - t0} seconds")


if __name__ == '__main__':
    br = BioReactor(compounds_path='data/compounds/drugs_paper_subset.csv',
                    output_path='results/results_example/',
                    organisms_path='data/organisms/organisms_to_use.tsv',
                    patterns_to_remove_path='data/patterns_to_remove/patterns.tsv',
                    molecules_to_remove_path='data/byproducts_to_remove/byproducts.tsv',
                    min_atom_count=5,
                    n_jobs=12)
    br.react()
