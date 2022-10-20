import logging

import click

from biocatalyzer import MSDataMatcher


@click.command()
@click.argument("ms_data",
                type=str,
                required=True,
                )
@click.argument("compounds_to_match",
                type=str,
                required=True,
                )
@click.argument("output_path",
                type=click.Path(),
                required=True,
                )
@click.option("--mode",
              "mode",
              type=click.Choice(['mass', 'mass_diff']),
              default='mass',
              show_default=True,
              help="The mode to use for the MS data matching (mass for ExactMass matching or mass_dif for ExactMass "
                   "differences matching).",
              )
@click.option("--tolerance",
              "tolerance",
              type=float,
              default=0.02,
              show_default=True,
              help="The mass tolerance to use when matching MS data.",
              )
def matcher_cli(ms_data,
                compounds_to_match,
                output_path,
                mode,
                tolerance):
    """Run the MSDataMatcher.

    Mandatory arguments:

        ms_data: Path to the file containing the MS data to use.

        compounds_to_match: Path to the file containing the compounds to match with the MS data.

        output_path: Path to the output directory.
    """
    ms = MSDataMatcher(ms_data_path=ms_data,
                       compounds_to_match=compounds_to_match,
                       output_path=output_path,
                       mode=mode,
                       tolerance=tolerance)
    logging.basicConfig(filename=f'{output_path}logging.log', level=logging.DEBUG)
    ms.generate_ms_results()


if __name__ == "__main__":
    pass
