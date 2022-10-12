# BioCatalyzer
...

## Resources:

Paper [Mapping human microbiome drug metabolism by gut bacteria and their genes](https://www.nature.com/articles/s41586-019-1291-3)

Supplementary data [here](https://www.nature.com/articles/s41586-019-1291-3#Sec53)!



## Command Line Interface

Usage of the script:

```bash
biocatalyzer <path_to_compounds> <output_path> [--reaction_rules=<FILE_PATH>] [--coreactants=<FILE_PATH>] [--patterns_to_remove=<FILE_PATH>] [--molecules_to_remove=<FILE_PATH>] [--min_atom_count=<INT>] [--n_jobs=<INT>]
```

| Argument            | Example           | Description                                                                                                                                                                                     | Default                                                                          |
|---------------------|-------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------|
| compounds           | file.tsv          | The path to the file containing the compounds to use as reactants.<sup>1</sup>                                                                                                                  |                                                                                  |
 | output_path         | output/directory/ | The path directory to save the results to.                                                                                                                                                      |                                                                                  |
| reaction_rules      | file.tsv          | The path to the file containing the reaction rules to use.<sup>2</sup>                                                                                                                          | [reaction_rules.csv](src/biocatalyzer/data/reactionrules/all_reaction_rules.tsv) |
| coreactants         | file.tsv          | The path to the file containing the coreactants to use.<sup>3</sup>                                                                                                                             | [coreactants.csv](src/biocatalyzer/data/coreactants/coreactants.tsv)             |
| organisms           | file.tsv          | The path to the file containing the organisms to use. Reaction Rules will be selected accordingly (select only rules associated with enzymes encoded by genes from this organisms).<sup>4</sup> | All reaction rules are used.                                                     |
| patterns_to_remove  | patterns.tsv      | The path to the file containing the patterns to remove from the products. <sup>5</sup>                                                                                                          | [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv)            |
| molecules_to_remove | molecules.tsv     | The path to the file containing the molecules to remove from the products. <sup>6</sup>                                                                                                         | [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv)      |
| min_atom_count      | 4                 | The minimum number of heavy atoms a product must have.                                                                                                                                          | 5                                                                                |
| n_jobs              | 6                 | The number of jobs to run in parallel (-1 uses all).                                                                                                                                            | 1                                                                                |

### Compounds

See [drugs.csv](src/biocatalyzer/data/compounds/drugs.csv)<sup>1</sup> for an example. 

The file must be tab separated and contain the following columns:
- `smiles` - the SMILES representation of the compounds;
- `compound_id` - the compounds identifiers.

### Reaction Rules

If you want to use your own Reaction Rules see [reaction_rules.csv](src/biocatalyzer/data/reactionrules/all_reaction_rules.tsv)<sup>2</sup> for an example.

The file must be tab separated and contain the following columns:

- `InternalID` - The ID of the Reaction Rule. # TODO: change the name of this column
- `Reactants` - The Reactants of the ReactionRule. Coreactants must be defined by their ID as in the Coreactants file.
The compound to match must be identifyed by the string 'Any'. The format must be: `coreactants1id;Any;coreactant2id`.
The order in which the reactants and the compound to match are defined is relevant and specific to the Reaction Rule.
If the Reaction Rules are mono-component (i.e. they do not contain any additional coreactant) the format must be: `Any`.
- `SMARTS` - The SMARTS representation of the Reaction Rule.
- `EC_Numbers` - The EC Numbers associated with the Reaction Rule.
- `Organisms` - The Organisms associated with the Reaction Rule.


### Coreactants

If you want to use your own Coreactants see [coreactants.csv](src/biocatalyzer/data/coreactants/all_coreactants.tsv)<sup>3</sup> for an example.

The file must be tab separated and contain the following columns:
- `smiles` - the SMILES representation of the compounds;
- `compound_id` - the compounds identifiers.


### Organisms

All organisms' identifiers are defined in: [https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html)

Example:

[hsa](https://www.genome.jp/kegg-bin/show_organism?org=hsa) is for *Homo sapiens* (human).

[eco](https://www.genome.jp/kegg-bin/show_organism?org=eco) is for *Escherichia coli K-12 MG1655*.

[sce](https://www.genome.jp/kegg-bin/show_organism?org=sce) is for *Saccharomyces cerevisiae (budding yeast)*.

<sup>4</sup> If you want to use your own organisms see [organisms.csv](src/biocatalyzer/data/organisms/organisms_to_use.csv) for an example.

The file must be tab separated and contain a column named `org_id` with the organisms' identifiers (KEGG identifiers).

### Patterns to remove

<sup>5</sup> If you want to use your own patterns to remove see [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv) for an example.

The file must be tab separated and contain a column named `smarts` with the SMARTS representation of the patterns to remove.

### Molecules to remove

<sup>6</sup> If you want to use your own molecules to remove see [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv) for an example.

The file must be tab separated and contain a column named `smiles` with the SMILES representation of the molecules to remove.
