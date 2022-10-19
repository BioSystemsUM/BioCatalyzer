import click


@click.command()
@click.argument("arg",
                type=str,
                required=True,
                )
def main(arg):
    print(arg)


if __name__ == "__main__":
    main()
