# BioReactor CLI

This CLI can be used to predicted metabolic products and associated enzymes for a set of input compounds.

**NOTE: if you want to use our full set of Reaction Rules download this [file](https://drive.google.com/file/d/1t2uYkKA8MjkIokSKNDESU27an1wW3CRK/view?usp=sharing) and provide its path in the `--reaction_rules` argument.**

## Command Line Interface

```bash
bioreactor_cli <PATH_TO_COMPOUNDS> <OUTPUT_DIRECTORY> [--neutralize=<BOOL>] [--reaction_rules=<FILE_PATH>] [--organisms=<FILE_PATH>] [--patterns_to_remove=<FILE_PATH>] [--molecules_to_remove=<FILE_PATH>] [--min_atom_count=<INT>] [--n_jobs=<INT>]
```

| Argument            | Example                                                 | Description                                                                                                                                                                                                                | Default                                                                                                                                                      |
|---------------------|---------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| compounds           | `file.tsv` or `"smile1;smiles2;smile3;etc"`             | The path to the file containing the compounds to use as reactants. Or ;-separated SMILES strings.                                                                                                                          |                                                                                                                                                              |
 | output_directory    | `output/directory/`                                     | The path directory to save the results to.                                                                                                                                                                                 |                                                                                                                                                              |
| neutralize          | `True` or `False`                                       | Whether to neutralize the compounds before predicting the products. In this case the new products will also be neutralized.                                                                                                | `False`                                                                                                                                                      |
| reaction_rules      | `file.tsv` or `None`                                    | The path to the file containing the reaction rules to use.                                                                                                                                                                 | [all_reaction_rules_forward_no_smarts_duplicates_sample.tsv](src/biocatalyzer/data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv) |
| organisms           | `file.tsv` or `"org_id1;org_id2;org_id3;etc"` or `None` | The path to the file containing the organisms to use. Or ;-separated organisms identifiers. Reaction Rules will be selected accordingly (select only rules associated with enzymes encoded by genes from this organisms).  | All reaction rules are used.                                                                                                                                 |
| patterns_to_remove  | `patterns.tsv` or `None`                                | The path to the file containing the patterns to remove from the products.                                                                                                                                                  | [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv)                                                                                        |
| molecules_to_remove | `molecules.tsv` or `None`                               | The path to the file containing the molecules to remove from the products.                                                                                                                                                 | [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv)                                                                                  |
| min_atom_count      | `4`                                                     | The minimum number of heavy atoms a product must have.                                                                                                                                                                     | `5`                                                                                                                                                          |
| n_jobs              | `6`                                                     | The number of jobs to run in parallel (-1 uses all).                                                                                                                                                                       | `1`                                                                                                                                                          |

More detailed information on these arguments and possible usages can be consulted in the main [README](README.md).

## Example:

```bash
bioreactor_cli compounds.tsv output_directory/ --neutralize=True --reaction_rules=reaction_rules.tsv --organisms=organisms.tsv --patterns_to_remove=patterns.tsv --molecules_to_remove=molecules.tsv --min_atom_count=4 --n_jobs=6
```