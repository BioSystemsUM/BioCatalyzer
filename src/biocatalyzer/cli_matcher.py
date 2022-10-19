import click


@click.command()
@click.argument("arg",
                type=str,
                required=True,
                )
def main(arg):
    print(f"{arg}_2")


if __name__ == "__main__":
    main()
