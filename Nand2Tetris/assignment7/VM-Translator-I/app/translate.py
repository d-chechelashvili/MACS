import os
from typing import List

import app.instruction_parser as parser
import app.translator as translator


def init_output_file(input_path: str) -> str:
    output_path = os.path.join(os.getcwd(), "output")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    _, file_name = os.path.split(input_path)
    output_name = file_name.split(".", 1)[0] + ".asm"
    output_path = os.path.join(output_path, output_name)
    with open(output_path, "w") as output_file:
        output_file.write("")
    return output_path


def write_in_file(instructions: List[str], input_path: str) -> None:
    output_path = init_output_file(input_path)
    with open(output_path, "a") as output_file:
        for instruction in instructions:
            output_file.write(instruction + "\n")


def translate(vm_file_name: str) -> None:
    instructions = parser.get_instructions(vm_file_name)
    binary_instructions = [
        translator.translate(instruction, vm_file_name) for instruction in instructions
    ]
    write_in_file(binary_instructions, vm_file_name)
