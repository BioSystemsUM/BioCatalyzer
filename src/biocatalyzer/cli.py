import click

from biocatalyzer.bioreactor import BioReactor
from biocatalyzer.matcher import MSDataMatcher


@click.command()
@click.argument("compounds",
                type=str,
                required=True,
                )
@click.argument("output_path",
                type=click.Path(),
                required=True,
                )
@click.option("--neutralize",
              "neutralize",
              type=bool,
              default=False,
              help="Whether to neutralize input compounds and newly generated compounds.")
@click.option("--organisms",
              "organisms",
              type=click.Path(exists=True),
              default=None,
              help="The path to the user defined file containing the organisms to filter the reaction rules.",
              )
@click.option("--patterns_to_remove",
              "patterns_to_remove",
              type=click.File('r'),
              default=None,
              show_default=True,
              help="A user defined file containing SMARTS patterns. Products that match a pattern will be removed.",
              )
@click.option("--molecules_to_remove",
              "molecules_to_remove",
              type=click.File('r'),
              default=None,
              show_default=True,
              help="A user defined file containing molecules encoded as SMILES to be removed from the products.",
              )
@click.option("--min_atom_count",
              "min_atom_count",
              type=int,
              default=5,
              show_default=True,
              help="The minimum atom count of a molecule (molecules with less atoms are removed from the products).",
              )
@click.option("--match_ms_data",
              "match_ms_data",
              type=bool,
              default=False,
              show_default=True,
              help="Whether to match the generated products with MS data.")
@click.option("--ms_data_path",
              "ms_data_path",
              type=click.File('r'),
              default=None,
              show_default=True,
              help="The path to the file containing the MS data to use.")
@click.option("--mode",
              "mode",
              type=click.Choice(['mass', 'mass_diff']),
              default='mass',
              show_default=True,
              help="The mode to use for the MS data matching (mass for ExactMass matching or mass_dif for ExactMass "
                   "differences matching).")
@click.option("--tolerance",
              "tolerance",
              type=float,
              default=0.02,
              show_default=True,
              help="The mass tolerance to use when matching MS data.",
              )
@click.option("--n_jobs",
              "n_jobs",
              type=int,
              default=1,
              show_default=True,
              help="The number of jobs to run in parallel.",
              )
def main(compounds,
         output_path,
         neutralize,
         organisms,
         patterns_to_remove,
         molecules_to_remove,
         min_atom_count,
         match_ms_data,
         ms_data_path,
         mode,
         tolerance,
         n_jobs):
    """Run the biocatalyzer.

    Mandatory arguments:

        compounds: Path to the file containing the compounds to use.

        output_path: Path to the output directory.
    """
    br = BioReactor(compounds_path=compounds,
                    output_path=output_path,
                    neutralize_compounds=neutralize,
                    organisms_path=organisms,
                    patterns_to_remove_path=patterns_to_remove,
                    molecules_to_remove_path=molecules_to_remove,
                    min_atom_count=min_atom_count,
                    n_jobs=n_jobs)
    br.react()

    if match_ms_data:
        if not ms_data_path:
            raise ValueError("The path to the MS data file is required when matching MS data.")

        ms = MSDataMatcher(ms_data_path=ms_data_path,
                           compounds_to_match=f"{output_path}/new_compounds.tsv",
                           output_path=output_path,
                           mode=mode,
                           tolerance=tolerance)

        ms.generate_ms_results()


if __name__ == "__main__":
    main()
