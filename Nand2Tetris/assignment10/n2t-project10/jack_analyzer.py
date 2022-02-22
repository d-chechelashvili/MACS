import click

from app import analyze


@click.command()
@click.option(
    "--target",
    required=True,
    help="The source Jack file or directory to be analyzed",
)
def cli(target: str) -> None:
    click.echo(f"Analyzing <{target}>")
    analyze(jack_file_or_directory_name=target)
    click.echo("Done!")


if __name__ == "__main__":
    cli()
