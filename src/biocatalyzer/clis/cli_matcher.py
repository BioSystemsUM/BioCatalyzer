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
def matcher_cli(ms_data,
                compounds_to_match,
                output_path,
                tolerance,
                n_jobs):
    """Run the MSDataMatcher.

    Mandatory arguments:

        ms_data: Path to the file containing the MS data to use.

        compounds_to_match: Path to the file containing the compounds to match with the MS data.

        output_path: Path to the output directory.
    """
    ms = MSDataMatcher(ms_data_path=ms_data,
                       compounds_to_match_path=compounds_to_match,
                       output_path=output_path,
                       tolerance=tolerance,
                       n_jobs=n_jobs)
    logging.basicConfig(filename=f'{output_path}_logging.log', level=logging.DEBUG)
    ms.generate_ms_results()


if __name__ == "__main__":
    pass
