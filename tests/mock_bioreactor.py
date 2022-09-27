import os
from typing import List

from biocatalyzer.bioreactor import BioReactor


class MockBioReactor(BioReactor):

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
