import click

from app import assemble


@click.command()
@click.option(
    "--asm-file",
    required=True,
    help="The source assembly file to be translated",
)
def cli(asm_file: str) -> None:
    click.echo(f"Hacking <{asm_file}>")
    assemble(asm_file)
    click.echo("Done!")


if __name__ == "__main__":
    cli()
