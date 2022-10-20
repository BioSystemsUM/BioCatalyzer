import click


@click.command()
@click.argument("arg",
                type=str,
                required=True,
                )
def matcher_cli(arg):
    print(f"{arg}_2")


if __name__ == "__main__":
    pass
