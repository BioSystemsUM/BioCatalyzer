from typing import Union

import pandas as pd

from biocatalyzer.chem import ChemUtils
from biocatalyzer.io_utils import Loaders


class MSDataMatcher:

    def __init__(self,
                 ms_data_path: str,
                 output_path: str,
                 mode: str = 'mass',
                 ms_field: str = 'Mass',
                 tolerance: float = 0.02,
                 reaction_rules: Union[pd.DataFrame, str] = None,
                 new_compounds: Union[pd.DataFrame, str] = None):
        self._ms_data = Loaders.load_ms_data(ms_data_path, ms_field)
        self._output_path = output_path
        self._mode = mode
        self._ms_field = ms_field
        self._tolerance = tolerance
        self._reaction_rules = reaction_rules
        self._new_compounds = new_compounds
        self._new_compounds['NewCompoundExactMass'] = \
            [ChemUtils.calc_exact_mass(m) for m in self._new_compounds.NewCompoundSMILES.values]

    def _match_masses(self):
        """
        Match the masses.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the matches.
        """
        pass

    def _match_mass_diff(self):
        pass

    def generate_ms_results(self):
        """
        Get the dataset with the matches.

        Returns
        -------
        pd.DataFrame:
            pandas dataframe with the matches.
        """
        if self._mode == 'mass':
            return self._match_masses()
        elif self._mode == 'mass_diff':
            return self._match_mass_diff()
        else:
            raise ValueError('The mode must be either "mass" or "mass_dif".')
