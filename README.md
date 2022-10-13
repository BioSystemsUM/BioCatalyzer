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
biocatalyzer <path_to_compounds> <output_path> [--reaction_rules=<FILE_PATH>] [--coreactants=<FILE_PATH>] [--patterns_to_remove=<FILE_PATH>] [--molecules_to_remove=<FILE_PATH>] [--min_atom_count=<INT>] [--masses=<FILE_PATH>or<str>] [--mass_tolerance=<float>] [--n_jobs=<INT>]
```

| Argument            | Example                                       | Description                                                                                                                                                                                                                           | Default                                                                          |
|---------------------|-----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------|
| compounds           | `file.tsv` or `"smile1;smiles2;smile3;etc"`   | The path to the file containing the compounds to use as reactants. Or ;-separated SMILES strings.<sup>1</sup>                                                                                                                         |                                                                                  |
 | output_path         | `output/directory/`                           | The path directory to save the results to.                                                                                                                                                                                            |                                                                                  |
| reaction_rules      | `file.tsv`                                    | The path to the file containing the reaction rules to use.<sup>2</sup>                                                                                                                                                                | [reaction_rules.csv](src/biocatalyzer/data/reactionrules/all_reaction_rules.tsv) |
| coreactants         | `file.tsv`                                    | The path to the file containing the coreactants to use.<sup>3</sup>                                                                                                                                                                   | [coreactants.csv](src/biocatalyzer/data/coreactants/all_coreactants.tsv)         |
| organisms           | `file.tsv` or `"org_id1;org_id2;org_id3;etc"` | The path to the file containing the organisms to use. Or ;-separated organisms identifiers. Reaction Rules will be selected accordingly (select only rules associated with enzymes encoded by genes from this organisms).<sup>4</sup> | All reaction rules are used.                                                     |
| patterns_to_remove  | `patterns.tsv`                                | The path to the file containing the patterns to remove from the products. <sup>5</sup>                                                                                                                                                | [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv)            |
| molecules_to_remove | `molecules.tsv`                               | The path to the file containing the molecules to remove from the products. <sup>6</sup>                                                                                                                                               | [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv)      |
| min_atom_count      | `4`                                           | The minimum number of heavy atoms a product must have.                                                                                                                                                                                | 5                                                                                |
| masses              | `masses.tsv` or `"mass1;mass2;mass2;etc"`     | The path to the file containing the masses of the compounds to use as reactants. Or ;-separated masses. <sup>7</sup>                                                                                                                  | None (Do not match any mass in specific.)                                        |
| mass_tolerance      | `0.02`                                        | The mass tolerance to use when matching masses.                                                                                                                                                                                       | 0.02                                                                             |
| n_jobs              | `6`                                           | The number of jobs to run in parallel (-1 uses all).                                                                                                                                                                                  | 1                                                                                |

### Compounds

See [drugs.csv](src/biocatalyzer/data/compounds/drugs.csv)<sup>1</sup> for an example. 

The file must be tab-separated and contain the following columns:
- `smiles` - the SMILES representation of the compounds;
- `compound_id` - the compounds identifiers.

Alternatively, the compounds can be passed as ;-separated string with the SMILES representations.

### Reaction Rules

If you want to use your own Reaction Rules see [reaction_rules.csv](src/biocatalyzer/data/reactionrules/all_reaction_rules.tsv)<sup>2</sup> for an example.

The file must be tab-separated and contain the following columns:

- `InternalID` - The ID of the Reaction Rule. # TODO: change the name of this column
- `Reactants` - The Reactants of the ReactionRule. Coreactants must be defined by their ID as in the Coreactants file.
The compound to match must be identifyed by the string 'Any'. The format must be: `coreactant1_id;Any;coreactant_id`.
The order in which the reactants and the compound to match are defined is relevant and specific to the Reaction Rule.
If the Reaction Rules are mono-component (i.e. they do not contain any additional coreactant) the format must be: `Any`.
- `SMARTS` - The SMARTS representation of the Reaction Rule.
- `EC_Numbers` - The EC Numbers associated with the Reaction Rule.
- `Organisms` - The Organisms associated with the Reaction Rule.


### Coreactants

If you want to use your own Coreactants see [coreactants.csv](src/biocatalyzer/data/coreactants/all_coreactants.tsv)<sup>3</sup> for an example.

The file must be tab-separated and contain the following columns:
- `smiles` - the SMILES representation of the compounds;
- `compound_id` - the compounds identifiers.


### Organisms

All organisms' identifiers are defined in: [https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html)

Example:

[hsa](https://www.genome.jp/kegg-bin/show_organism?org=hsa) is for *Homo sapiens* (human).

[eco](https://www.genome.jp/kegg-bin/show_organism?org=eco) is for *Escherichia coli K-12 MG1655*.

[sce](https://www.genome.jp/kegg-bin/show_organism?org=sce) is for *Saccharomyces cerevisiae (budding yeast)*.

<sup>4</sup> If you want to use your own organisms see [organisms.csv](src/biocatalyzer/data/organisms/organisms_to_use.tsv) for an example.

The file must be tab-separated and contain a column named `org_id` with the organisms' identifiers (KEGG identifiers).

Alternatively, the organisms can be passed as ;-separated string with the organisms identifiers.

### Patterns to remove

<sup>5</sup> If you want to use your own patterns to remove see [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv) for an example.

The file must be tab-separated and contain a column named `smarts` with the SMARTS representation of the patterns to remove.

### Molecules to remove

<sup>6</sup> If you want to use your own molecules to remove see [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv) for an example.

The file must be tab-separated and contain a column named `smiles` with the SMILES representation of the molecules to remove.

### Masses

You can provide a tab-separated file containing the masses to match the new products. See [masses.tsv](src/biocatalyzer/data/masses_to_match/masses.tsv) for an example.

The file must be tab-separated and contain a column named `mass` with the masses to match.

Alternatively, the masses can be passed as ;-separated string with the masses to match.

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