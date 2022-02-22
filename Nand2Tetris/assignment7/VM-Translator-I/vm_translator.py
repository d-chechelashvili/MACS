import click

from app import translate


@click.command()
@click.option(
    "--vm-file",
    required=True,
    help="The source VM file to be translated",
)
def cli(vm_file: str) -> None:
    click.echo(f"Translating <{vm_file}>")
    translate(vm_file)
    click.echo("Done!")


if __name__ == "__main__":
    cli()
