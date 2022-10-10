# BioCatalyzer
...

## Command Line Interface

Usage of the script:

```bash
python cli.py <path_to_compounds> <output_path> [--reaction_rules=<FILE_PATH>] [--coreactants=<FILE_PATH>] [--patterns_to_remove=<FILE_PATH>] [--molecules_to_remove=<FILE_PATH>] [--min_atom_count=<INT>] [--n_jobs=<INT>]
```

| Argument            | Example            | Description                                                                             | Default                                   |
|---------------------|--------------------|-----------------------------------------------------------------------------------------|-------------------------------------------|
| compounds           | file.csv           | The path to the file containing the compounds to use as reactants.<sup>1</sup>          |                                           |
 | output_path         | output/directory/  | The path directory to save the results to.                                              |                                           |
| patterns_to_remove  | patterns.tsv       | The path to the file containing the patterns to remove from the products. <sup>4</sup>  | data/patterns_to_remove/patterns.tsv      |
| molecules_to_remove | molecules.tsv      | The path to the file containing the molecules to remove from the products. <sup>5</sup> | data/byproducts_to_remove/byproducts.tsv  |
| min_atom_count      | 4                  | The minimum number of heavy atoms a product must have.                                  | 5                                         |
| organisms           | hsa                | Organisms to use (select ony enzymes encoded by genes from this organisms).             | all                                       |
| n_jobs              | 12                 | The number of jobs to run in parallel.                                                  | 1                                         |

<sup>1</sup> See [drugs.csv](src/biocatalyzer/data/compounds/drugs.csv) for an example.

<sup>4</sup> See [patterns.tsv](src/biocatalyzer/data/patterns_to_remove/patterns.tsv) for an example.

<sup>5</sup> See [byproducts.tsv](src/biocatalyzer/data/byproducts_to_remove/byproducts.tsv) for an example.

## Organisms:

All organisms' ids are defined in: [https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html)

Example:

[hsa](https://www.genome.jp/kegg-bin/show_organism?org=hsa) is for *Homo sapiens* (human).

[eco](https://www.genome.jp/kegg-bin/show_organism?org=eco) is for *Escherichia coli K-12 MG1655*.

[sce](https://www.genome.jp/kegg-bin/show_organism?org=sce) is for *Saccharomyces cerevisiae (budding yeast)*.

etc...