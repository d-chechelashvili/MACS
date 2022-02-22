import click

from app import compile


@click.command()
@click.option(
    "--target",
    required=True,
    help="The source Jack file or directory to be compiled",
)
def cli(target: str) -> None:
    click.echo(f"Compiling <{target}>")
    compile(jack_file_or_directory_name=target)
    click.echo("Done!")


if __name__ == "__main__":
    cli()
