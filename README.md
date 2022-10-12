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
| organisms           | file.tsv          | The path to the file containing the organisms to use. Reaction Rules will be selected accordingly (select only rules associated with enzymes encoded by genes from this organisms).<sup>3</sup> | All reaction rules are used.                                                     |
| patterns_to_remove  | patterns.tsv      | The path to the file containing the patterns to remove from the products. <sup>3</sup>                                                                                                          | [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv)            |
| molecules_to_remove | molecules.tsv     | The path to the file containing the molecules to remove from the products. <sup>4</sup>                                                                                                         | [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv)      |
| min_atom_count      | 4                 | The minimum number of heavy atoms a product must have.                                                                                                                                          | 5                                                                                |
| n_jobs              | 6                 | The number of jobs to run in parallel (-1 uses all).                                                                                                                                            | 1                                                                                |

<sup>1</sup> See [drugs.csv](src/biocatalyzer/data/compounds/drugs.csv) for an example.

<sup>2</sup> See [reaction_rules.csv](src/biocatalyzer/data/reactionrules/all_reaction_rules.tsv) for an example.

<sup>3</sup> See [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv) for an example.

<sup>4</sup> See [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv) for an example.

## Organisms:

All organisms' ids are defined in: [https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html)

Example:

[hsa](https://www.genome.jp/kegg-bin/show_organism?org=hsa) is for *Homo sapiens* (human).

[eco](https://www.genome.jp/kegg-bin/show_organism?org=eco) is for *Escherichia coli K-12 MG1655*.

[sce](https://www.genome.jp/kegg-bin/show_organism?org=sce) is for *Saccharomyces cerevisiae (budding yeast)*.

etc...