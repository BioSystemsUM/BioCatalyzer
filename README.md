# BioCatalyzer

BioCatalyzer is a python tool that predicts enzymatic metabolism products using a rule-based approach.

BioCatalyzer is implemented as a Command Line Interface that takes as input a set of compounds represented as SMILES 
strings and outputs a set of predicted metabolic products and associated enzymes.

<!---
## Resources:

Paper [Mapping human microbiome drug metabolism by gut bacteria and their genes](https://www.nature.com/articles/s41586-019-1291-3)

Supplementary data [here](https://www.nature.com/articles/s41586-019-1291-3#Sec53)!
--->

## Instalation (Not working yet!)

Installing from Pypi package repository:

`pip install biocatalyzer`

Installing from GitHub:

1. clone the repository: `git clone https://github.com/jcorreia11/BioCatalyzer.git`

2. run: `python setup.py install`

## Command Line Interface

```bash
biocatalyzer <PATH_TO_COMPOUNDS> <OUTPUT_DIRECTORY> [--neutralize=<BOOL>] [--reaction_rules=<FILE_PATH>] [--organisms=<FILE_PATH>] [--patterns_to_remove=<FILE_PATH>] [--molecules_to_remove=<FILE_PATH>] [--min_atom_count=<INT>] [--match_ms_data=<BOOL>] [--ms_data_path=<FILE_PATH>] [--mode=<STR>] [--tolerance=<FLOAT>] [--n_jobs=<INT>]
```

| Argument            | Example                                                 | Description                                                                                                                                                                                                                           | Default                                                                                                                                                      |
|---------------------|---------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| compounds           | `file.tsv` or `"smile1;smiles2;smile3;etc"`             | The path to the file containing the compounds to use as reactants. Or ;-separated SMILES strings.<sup>1</sup>                                                                                                                         |                                                                                                                                                              |
 | output_directory    | `output/directory/`                                     | The path directory to save the results to.                                                                                                                                                                                            |                                                                                                                                                              |
| neutralize          | `True` or `False`                                       | Whether to neutralize the compounds before predicting the products. In this case the new products will also be neutralized.                                                                                                           | `False`                                                                                                                                                      |
| reaction_rules      | `file.tsv` or `None`                                    | The path to the file containing the reaction rules to use.<sup>2</sup>                                                                                                                                                                | [all_reaction_rules_forward_no_smarts_duplicates_sample.tsv](src/biocatalyzer/data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv) |
| organisms           | `file.tsv` or `"org_id1;org_id2;org_id3;etc"` or `None` | The path to the file containing the organisms to use. Or ;-separated organisms identifiers. Reaction Rules will be selected accordingly (select only rules associated with enzymes encoded by genes from this organisms).<sup>3</sup> | All reaction rules are used.                                                                                                                                 |
| patterns_to_remove  | `patterns.tsv` or `None`                                | The path to the file containing the patterns to remove from the products. <sup>4</sup>                                                                                                                                                | [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv)                                                                                        |
| molecules_to_remove | `molecules.tsv` or `None`                               | The path to the file containing the molecules to remove from the products. <sup>5</sup>                                                                                                                                               | [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv)                                                                                  |
| min_atom_count      | `4`                                                     | The minimum number of heavy atoms a product must have.                                                                                                                                                                                | `5`                                                                                                                                                          |
| match_ms_data       | `True` or `False`                                       | Whether to match the predicted products to the MS data.                                                                                                                                                                               | `False`                                                                                                                                                      |
| ms_data_path        | `ms_data.tsv`                                           | The path to the file containing the MS data. <sup>6</sup>                                                                                                                                                                             | `None`                                                                                                                                                       |
| mode                | `mass` or `mass_diff`                                   | The mode to use when matching the predicted products to the MS data.                                                                                                                                                                  | `mass`                                                                                                                                                       |
| tolerance           | `0.02`                                                  | The mass tolerance to use when matching masses.                                                                                                                                                                                       | `0.02`                                                                                                                                                       |
| n_jobs              | `6`                                                     | The number of jobs to run in parallel (-1 uses all).                                                                                                                                                                                  | `1`                                                                                                                                                          |

### Compounds

See [drugs.csv](src/biocatalyzer/data/compounds/drugs.csv)<sup>1</sup> for an example. 

The file must be tab-separated and contain the following columns:
- `smiles` - the SMILES representation of the compounds;
- `compound_id` - the compounds identifiers.

Alternatively, the compounds can be passed as ;-separated string with the SMILES representations.

### Output directory

The output path must be a directory. The results will be saved in the following files:
- `new_compouds.tsv` - the predicted products;
- `matches.tsv` (if `match_ms_data` is set to `True`) - the matches between the predicted products and the MS data;

### Neutralize

If set to `True`, the compounds will be neutralized before predicting the products. In this case the new products will 
also be neutralized.

### Reaction Rules

See [all_reaction_rules_forward_no_smarts_duplicates_sample.tsv](src/biocatalyzer/data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv)<sup>2</sup> 
for an example.

The file must be tab-separated and contain the following columns:

- `InternalID` - The ID of the Reaction Rule. # TODO: change the name of this column
- `Reactants` - The Reactants of the ReactionRule. Coreactants must be defined by their ID as in the Coreactants file.
The compound to match must be identifyed by the string 'Any'. The format must be: `coreactant1_id;Any;coreactant_id`.
The order in which the reactants and the compound to match are defined is relevant and specific to the Reaction Rule.
If the Reaction Rules are mono-component (i.e. they do not contain any additional coreactant) the format must be: `Any`.
- `SMARTS` - The SMARTS representation of the Reaction Rule.
- `EC_Numbers` - The EC Numbers associated with the Reaction Rule.
- `Organisms` - The Organisms associated with the Reaction Rule.

To use our complete set of Reaction Rules please download the following 
[file](https://drive.google.com/file/d/1t2uYkKA8MjkIokSKNDESU27an1wW3CRK/view?usp=sharing).

You can directly use this file by providing the path to it as the value of the `reaction_rules` parameter.

### Organisms

All organisms' identifiers are defined in: 
[https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html) are allowed.

Example:

[hsa](https://www.genome.jp/kegg-bin/show_organism?org=hsa) is for *Homo sapiens* (human).

[eco](https://www.genome.jp/kegg-bin/show_organism?org=eco) is for *Escherichia coli K-12 MG1655*.

[sce](https://www.genome.jp/kegg-bin/show_organism?org=sce) is for *Saccharomyces cerevisiae (budding yeast)*.

If you want to use your own organisms see 
[organisms.csv](src/biocatalyzer/data/organisms/organisms_to_use.tsv)<sup>3</sup> for an example.

The file must be tab-separated and contain a column named `org_id` with the organisms' identifiers (KEGG identifiers).

Alternatively, the organisms can be passed as ;-separated string with the organisms identifiers.

### Patterns to remove

If you want to use your own patterns to remove see 
[patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv)<sup>4</sup> for an example.

The file must be tab-separated and contain a column named `smarts` with the SMARTS representation of the patterns to remove.

### Molecules to remove

If you want to use your own molecules to remove see 
[byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv)<sup>5</sup> for an example.

The file must be tab-separated and contain a column named `smiles` with the SMILES representation of the molecules to remove.

### Match MS data

If set to `True`, the predicted products will be matched to the MS data.

In this case the `ms_data_path` must be set.

### MS data path

See [ms_data.tsv](src/biocatalyzer/data/ms_data_example/ms_data_paper.tsv)<sup>6</sup> for an example.

The file must be tab-separated and contain the following columns:
- `ParentCompound` - the parent/original compound identifiers.
- `ParentCompoundSmiles` - the SMILES representation of the compounds (optional).
- `Mass` or `MassDiff` - depending on the selected `mode`, the mass or mass difference of the molecule.

### Mode

The mode to use when matching the predicted products to the MS data.

If set to `mass`, the `Mass` column will be used. This will match the predicted products exact mass to the MS 
data provided in the `Mass` column.

If set to `mass_diff`, the `MassDiff` column will be used. This will match the predicted products mass difference to 
the ParentDrug as provided in the MS data `MassDiff` column.

### Mass Tolerance

The mass tolerance (`float`) to use when matching masses. Masses between `mass - mass_tolerance` and `mass + mass_tolerance` will be considered as a match.

### Number of jobs

The number of jobs to run in parallel. If `-1` is passed, all available cores will be used.

## Cite

Manuscript under preparation!

### Credits and License

Developed at Centre of Biological Engineering, University of Minho and EMBL ...

This project has received funding from the Portuguese FCT and EMBL (?).

Released under an MIT License. # TODO: check if licence is in accordance with packages/data used