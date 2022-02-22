import click

from app import translate


@click.command()
@click.option(
    "--target",
    required=True,
    help="The source VM file or directory to be translated",
)
def cli(target: str) -> None:
    click.echo(f"Translating <{target}>")
    translate(vm_file_or_directory_name=target)
    click.echo("Done!")


if __name__ == "__main__":
    cli()
