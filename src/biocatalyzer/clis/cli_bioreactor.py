import logging
import os

import click

from biocatalyzer import BioReactor

DATA_FILES = os.path.dirname(__file__)


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
              help="Whether to neutralize input compounds and newly generated compounds.",
              )
@click.option("--reaction_rules",
              "reaction_rules",
              type=str,
              default='default',
              show_default=True,
              help="Path to reaction rules file.",
              )
@click.option("--organisms",
              "organisms",
              type=str,
              default=None,
              help="The path to the user defined file containing the organisms to filter the reaction rules.",
              )
@click.option("--patterns_to_remove",
              "patterns_to_remove",
              type=str,
              default='default',
              show_default=True,
              help="A user defined file containing SMARTS patterns. Products that match a pattern will be removed.",
              )
@click.option("--molecules_to_remove",
              "molecules_to_remove",
              type=str,
              default='default',
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
@click.option("--n_jobs",
              "n_jobs",
              type=int,
              default=1,
              show_default=True,
              help="The number of jobs to run in parallel.",
              )
def bioreactor_cli(compounds,
                   output_path,
                   neutralize,
                   reaction_rules,
                   organisms,
                   patterns_to_remove,
                   molecules_to_remove,
                   min_atom_count,
                   n_jobs):
    """Run the BioCatalyzer.

    Mandatory arguments:

        compounds: Path to the file containing the compounds to use.

        output_path: Path to the output directory.
    """
    if reaction_rules is None:
        reaction_rules = os.path.join(
            DATA_FILES, '../data/reactionrules/all_reaction_rules_forward_no_smarts_duplicates_sample.tsv')
    br = BioReactor(compounds_path=compounds,
                    output_path=output_path,
                    reaction_rules_path=reaction_rules,
                    neutralize_compounds=neutralize,
                    organisms_path=organisms,
                    patterns_to_remove_path=patterns_to_remove,
                    molecules_to_remove_path=molecules_to_remove,
                    min_atom_count=min_atom_count,
                    n_jobs=n_jobs)
    logging.basicConfig(filename=f'{output_path}logging.log', level=logging.DEBUG)
    br.react()


if __name__ == "__main__":
    pass
