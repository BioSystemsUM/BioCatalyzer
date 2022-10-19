from typing import List


def _correct_number_of_parenthesis(smiles_list: List[str]):
    """
    Corrects the number of parenthesis in a list of SMILES strings.
    Sometimes the react method returns a SMILES string with an incorrect number of parenthesis.
    This method corrects that issue.

    Parameters
    ----------
    smiles_list: List[str]
        The list of SMILES strings to correct the parenthesis.
    """
    corrected_smiles = []
    for p in smiles_list:
        # deal with cases where invalid number of parentheses are generated
        if (p.count('(') + p.count(')')) % 2 != 0:
            if p[0] == '(':
                p = p[1:]
            elif p[-1] == ')':
                p = p[:-1]
        corrected_smiles.append(p)
    return corrected_smiles
