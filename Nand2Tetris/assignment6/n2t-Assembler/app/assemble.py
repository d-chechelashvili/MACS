import os
from typing import List

from app.instructionParser import get_instructions
from app.SymbolTable import SymbolTable
from app.translator import translate


def init_output_file(input_path: str) -> str:
    output_path = os.path.join(os.getcwd(), "output")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    _, file_name = os.path.split(input_path)
    output_name = file_name.split(".", 1)[0] + ".hack"
    output_path = os.path.join(output_path, output_name)
    with open(output_path, "w") as output_file:
        output_file.write("")
    return output_path


def write_in_file(instructions: List[str], input_path: str) -> None:
    output_path = init_output_file(input_path)
    with open(output_path, "a") as output_file:
        for instruction in instructions:
            output_file.write(instruction + "\n")


def assemble(asm_file: str) -> None:
    symbol_table = SymbolTable()
    instructions = get_instructions(asm_file, symbol_table)
    binary_instructions = [
        translate(instruction, symbol_table) for instruction in instructions
    ]
    write_in_file(binary_instructions, asm_file)
