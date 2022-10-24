def _correct_number_of_parenthesis(smiles: str):
    """
    Corrects the number of parenthesis of a SMILES string.
    Sometimes the react method returns a SMILES string with an incorrect number of parenthesis.
    This method corrects that issue.

    Parameters
    ----------
    smiles: str
        The SMILES string to correct the parenthesis.

    Returns
    -------
    str
        The corrected SMILES string.
    """
    # deal with cases where invalid number of parentheses are generated
    if (smiles.count('(') + smiles.count(')')) % 2 != 0:
        if smiles[0] == '(':
            smiles = smiles[1:]
        elif smiles[-1] == ')':
            smiles = smiles[:-1]
    return smiles
