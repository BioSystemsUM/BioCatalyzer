# MSDataMatcher CLI

This CLI can be used to match the predicted metabolic products with experimental MS data.

## Command Line Interface

```bash
matcher_cli <MS_DATA> <COMPOUNDS_TO_MATCH> <OUTPUT_DIRECTORY> [--mode=<STR>] [--tolerance=<FLOAT>]
```

| Argument                                 | Example               | Description                                                                        | Default                                                                                                                                            |
|------------------------------------------|-----------------------|------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| ms_data <MS_DATA>                        | `ms_data.tsv`         | The path to the file containing the MS data.                                       |                                                                                                                                                    |
| compounds_to_match <COMPOUNDS_TO_MATCH>  | `file.tsv`            | The path to the file containing the predicted compounds to match with the MS data. |                                                                                                                                                    |
 | output_directory <OUTPUT_DIRECTORY>      | `output/directory/`   | The path directory to save the results to.                                         |                                                                                                                                                    |
| mode                                     | `mass` or `mass_diff` | The mode to use when matching the predicted products to the MS data.               | `mass`                                                                                                                                             |
| tolerance                                | `0.02`                | The mass tolerance to use when matching masses.                                    | `0.02`                                                                                                                                             |

More detailed information on these arguments and possible usages can be consulted in the main [README](README.md).

## Example:

```bash
matcher_cli ms_data.tsv compounds_to_match.tsv output_directory/ --mode=mass --tolerance=0.0015
```
