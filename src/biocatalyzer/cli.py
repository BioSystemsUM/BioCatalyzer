import click

from biocatalyzer.bioreactor import BioReactor


@click.command()
@click.argument("compounds",
                type=click.Path(exists=True),
                required=True,
                )
@click.argument("output_path",
                type=click.Path(),
                required=True,
                )
@click.option("--reaction_rules",
              "reaction_rules",
              type=click.Path(exists=True),
              default="data/reactionrules/all_reaction_rules.tsv",
              show_default=True,
              help="Path to the file containing the reaction rules to use.",
              )
@click.option("--coreactants",
              "coreactants",
              type=click.Path(exists=True),
              default="data/coreactants/all_coreactants.tsv",
              show_default=True,
              help="Path to the file containing the coreactants to use.",
              )
@click.option("--patterns_to_remove",
              "patterns_to_remove",
              type=click.Path(exists=True),
              default="data/patterns_to_remove/patterns.tsv",
              show_default=True,
              help="A file containing SMARTS patterns. Products that match a pattern will be removed.",
              )
@click.option("--molecules_to_remove",
              "molecules_to_remove",
              type=click.Path(exists=True),
              default="data/byproducts_to_remove/byproducts.tsv",
              show_default=True,
              help="A file containing molecules encoded as SMILES to be removed from the products.",
              )
@click.option("--min_atom_count",
              "min_atom_count",
              type=int,
              default=4,
              show_default=True,
              help="The minimum atom count of a molecule (molecules with less atoms are removed from the products).",
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
         reaction_rules,
         coreactants,
         patterns_to_remove,
         molecules_to_remove,
         min_atom_count,
         n_jobs):
    """Run the biocatalyzer.

    Mandatory arguments:
        compounds: Path to the file containing the compounds to use.
        output_path: Path to the output directory.
    """
    br = BioReactor(compounds_path=compounds,
                    output_path=output_path,
                    reaction_rules_path=reaction_rules,
                    coreactants_path=coreactants,
                    patterns_to_remove_path=patterns_to_remove,
                    molecules_to_remove_path=molecules_to_remove,
                    min_atom_count=min_atom_count,
                    n_jobs=n_jobs)
    br.react()


if __name__ == "__main__":
    main()