import multiprocessing
import os
import time
from typing import List

import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles

from biocatalyzer._utils import ChemUtils
from biocatalyzer.io_utils import Loaders


class BioReactor:
    """
    Main class of the BioCatalyzer package.
    Performs the transformation of the reactants into products using reaction rules.
    """

    def __init__(self,
                 compounds_path: str,
                 reaction_rules_path: str,
                 coreactants_path: str,
                 molecules_to_remove_path: str,
                 patterns_to_remove_path: str,
                 min_atom_count: int,
                 output_path: str,
                 n_jobs: int = 1):
        """
        Initialize the BioReactor class.

        Parameters
        ----------
        compounds_path: str
            The path to the file containing the compounds to use as reactants.
        reaction_rules_path: str
            The path to the file containing the reaction rules to use.
        coreactants_path: str
            The path to the file containing the coreactants to use.
        molecules_to_remove_path: str
            The path to the file containing the molecules to remove from the products.
        patterns_to_remove_path: str
            The path to the file containing the patterns to remove from the products.
        output_path: str
            The path to the output directory.
        n_jobs: int
            The number of jobs to run in parallel.
        """
        # silence RDKit logger
        RDLogger.DisableLog('rdApp.*')
        self._verify_files([compounds_path, reaction_rules_path, coreactants_path])
        self._set_output_path(output_path)
        self._compounds = Loaders.load_compounds(compounds_path)
        self._reaction_rules = Loaders.load_reaction_rules(reaction_rules_path)
        self._coreactants = Loaders.load_coreactants(coreactants_path)
        self._molecules_to_remove = Loaders.load_byproducts_to_remove(molecules_to_remove_path)
        self._patterns_to_remove = Loaders.load_patterns_to_remove(patterns_to_remove_path)
        self._min_atom_count = min_atom_count
        self._output_path = output_path
        self._n_jobs = n_jobs
        self._new_products = set()

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
        reactants_ids = self._reaction_rules[self._reaction_rules.smarts == smarts].coreactants_ids.values[0]
        reactants = self._set_reactants(reactants_ids, smiles)
        results = ChemUtils.react(reactants, smarts)
        if len(results) > 0:
            smiles_id = self._compounds[self._compounds.smiles == smiles].compound_id.values[0]
            smarts_id = self._reaction_rules[self._reaction_rules.smarts == smarts].rule_id.values[0]
            reaction_info = self._reaction_rules[self._reaction_rules.smarts == smarts].reaction_info.values[0]
            with open(self._output_path + '/results.tsv', 'a') as rs, \
                    open(self._output_path + '/new_compounds.tsv', 'a') as nc:
                i = 0
                for result in results:
                    rs.write(f"{smiles_id}\t{smarts_id}\t{result}\t{reaction_info}\n")
                    products = result.split('>')[-1].split('.')
                    for p in products:
                        # deal with cases where invalid number of parentheses are present
                        if (p.count('(') + p.count(')')) % 2 != 0:
                            if p[0] == '(':
                                p = p[1:]
                            elif p[-1] == ')':
                                p = p[:-1]
                        if p not in self._new_products:
                            # because multiprocessing is used, this does not guarantee that the same product will not
                            # be added multiple times
                            self._new_products.add(p)
                            if not self._match_byproducts(p) and not self._match_patterns(p) \
                                    and self._min_atom_count_filter(p):
                                nc.write(f"{smiles_id}_{smarts_id}_{i}\t{p}\n")
                                i += 1

    def react(self):
        """
        Transform reactants into products using the reaction rules.
        """
        t0 = time.time()
        with open(self._output_path + '/results.tsv', 'w') as rs, \
                open(self._output_path + '/new_compounds.tsv', 'w') as nc:
            rs.write('id_mol\tid_rule\tnew_reaction_smiles\treaction_info\n')
            nc.write('new_compound_id\tnew_compound_smiles\n')
        for compound in self._compounds.smiles:
            with multiprocessing.Pool(self._n_jobs) as pool:
                pool.starmap(self._react_single, zip([compound]*self._reaction_rules.shape[0],
                                                     self._reaction_rules.smarts))
        # TODO: change this later
        results = pd.read_csv(self._output_path + '/new_compounds.tsv', sep='\t')
        results.drop_duplicates(inplace=True, subset=['new_compound_smiles'])
        results.to_csv(self._output_path + '/new_compounds.tsv', sep='\t', index=False)

        t1 = time.time()
        print(f"{results.shape[0]} unique new compounds generated!")
        print(f"Time elapsed: {t1 - t0} seconds")


if __name__ == '__main__':
    br = BioReactor(compounds_path='data/compounds/drugs.csv',
                    reaction_rules_path='data/reactionrules/all_reaction_rules.tsv',
                    coreactants_path='data/coreactants/all_coreactants.tsv',
                    patterns_to_remove_path='data/patterns_to_remove/patterns.tsv',
                    molecules_to_remove_path='data/byproducts_to_remove/byproducts.tsv',
                    min_atom_count=5,
                    output_path='results',
                    n_jobs=12)
    br.react()
