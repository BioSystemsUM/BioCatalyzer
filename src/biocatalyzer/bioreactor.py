import multiprocessing
import os
import time
from typing import List

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
                 reaction_rules_path: str = None,
                 coreactants_path: str = None,
                 organisms_path: str = None,
                 molecules_to_remove_path: str = None,
                 patterns_to_remove_path: str = None,
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
        reaction_rules_path: str
            The path to the file containing the reaction rules to use.
        coreactants_path: str
            The path to the file containing the coreactants to use.
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
        self._reaction_rules_path = reaction_rules_path
        self._coreactants_path = coreactants_path
        self._organisms_path = organisms_path
        self._molecules_to_remove_path = molecules_to_remove_path
        self._patterns_to_remove_path = patterns_to_remove_path
        self._set_up_files()
        self._orgs = Loaders.load_organisms(self._organisms_path)
        self._reaction_rules = Loaders.load_reaction_rules(self._reaction_rules_path, orgs=self._orgs)
        self._set_output_path(self._output_path)
        self._compounds = Loaders.load_compounds(self._compounds_path)
        self._coreactants = Loaders.load_coreactants(self._coreactants_path)
        self._molecules_to_remove = Loaders.load_byproducts_to_remove(self._molecules_to_remove_path)
        self._patterns_to_remove = Loaders.load_patterns_to_remove(self._patterns_to_remove_path)
        self._min_atom_count = min_atom_count
        if n_jobs == -1:
            self._n_jobs = multiprocessing.cpu_count()
        else:
            self._n_jobs = n_jobs

    def _set_up_files(self):
        self.DATA_FILES = os.path.join(os.path.dirname(__file__), 'data')
        if not self._reaction_rules_path:
            self._reaction_rules_path = os.path.join(self.DATA_FILES, 'reactionrules/all_reaction_rules_forward_no_smarts_duplicates.tsv')
        if not self._coreactants_path:
            self._coreactants_path = os.path.join(self.DATA_FILES, 'coreactants/all_coreactants.tsv')
        if not self._molecules_to_remove_path:
            self._molecules_to_remove_path = os.path.join(self.DATA_FILES, 'byproducts_to_remove/byproducts.tsv')
        if not self._patterns_to_remove_path:
            self._patterns_to_remove_path = os.path.join(self.DATA_FILES, 'patterns_to_remove/patterns.tsv')
        if self._organisms_path:
            self._verify_files([self._compounds_path, self._reaction_rules_path, self._coreactants_path,
                                self._organisms_path, self._molecules_to_remove_path, self._patterns_to_remove_path])
        else:
            self._verify_files([self._compounds_path, self._reaction_rules_path, self._coreactants_path,
                                self._molecules_to_remove_path, self._patterns_to_remove_path])

    @staticmethod
    def _verify_files(paths: List[str]):
        """
        Verify that the provided paths to the files exist.

        Parameters
        ----------
        paths: List[str]
            The paths to the files to verify.
        """
        for path in paths:
            if path:
                if not os.path.exists(path):
                    raise FileNotFoundError(f"File {path} not found.")

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

    def _set_reactants(self, reactants: str, compound: str):
        """
        Gets the coreactants in the correct order to be processed by the reaction rule.

        Parameters
        ----------
        reactants: str
            The coreactants ids to use.
        compound: str
            The main compound to use as reactant.

        Returns
        -------
        List[Mol]
            The coreactants in the correct order as Mol objects.
        """
        reactants_list = []
        reactants = reactants.split(';')
        if len(reactants) > 1:
            for r in reactants:
                if r == 'Any':
                    reactants_list.append(compound)
                else:
                    reactants_list.append(self._coreactants[self._coreactants.compound_id == r].smiles.values[0])
        else:
            if reactants[0] == 'Any':
                reactants_list.append(compound)
        return reactants_list

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
        nc = pd.DataFrame(columns=['ResultID', 'NewCompoundID', 'NewCompoundSMILES', 'EC_Numbers'])
        reactants_ids = self._reaction_rules[self._reaction_rules.SMARTS == smarts].Reactants.values[0]
        reactants = self._set_reactants(reactants_ids, smiles)
        results = ChemUtils.react(reactants, smarts)
        if len(results) > 0:
            smiles_id = self._compounds[self._compounds.smiles == smiles].compound_id.values[0]
            smarts_id = self._reaction_rules[self._reaction_rules.SMARTS == smarts].InternalID.values[0]
            with open(self._output_path + '/results.tsv', 'a') as rs:
                i = 0
                for i, result in enumerate(results):
                    id_result = f"{smiles_id}_{smarts_id}_{i}"
                    rs.write(f"{id_result}\t{smiles_id}\t{smarts_id}\t{result}\n")
                    products = result.split('>')[-1].split('.')
                    for p in products:
                        # deal with cases where invalid number of parentheses are generated
                        if (p.count('(') + p.count(')')) % 2 != 0:
                            if p[0] == '(':
                                p = p[1:]
                            elif p[-1] == ')':
                                p = p[:-1]
                        if p not in nc.NewCompoundSMILES.values:
                            if not self._match_byproducts(p) and not self._match_patterns(p) \
                                    and self._min_atom_count_filter(p):
                                ecs = self._get_ec_numbers(smarts_id)
                                nc.loc[len(nc)] = [id_result, f"{id_result}_{i}", p, ecs]
                                i += 1
        return nc

    def react(self):
        """
        Transform reactants into products using the reaction rules.
        """
        t0 = time.time()
        with open(self._output_path + '/results.tsv', 'w') as rs:
            rs.write('ResultID\tCompoundID\tRuleID\tNewReactionSmiles\n')
        for compound in self._compounds.smiles:
            with multiprocessing.Pool(self._n_jobs) as pool:
                results = pool.starmap(self._react_single, zip([compound]*self._reaction_rules.shape[0],
                                                               self._reaction_rules.SMARTS))

        results = pd.concat(results)
        results = results.drop_duplicates(subset=['NewCompoundSMILES'])
        results.to_csv(self._output_path + '/new_compounds.tsv', sep='\t', index=False)
        t1 = time.time()
        print(f"{results.shape[0]} unique new compounds generated!")
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
