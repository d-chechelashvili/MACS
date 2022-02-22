import os
from typing import List

import app.instruction_parser as parser
import app.translator as translator
from app.translator_II import get_boot_code


def init_output_file(input_path: str) -> str:
    output_path = os.path.join(os.getcwd(), "output")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    _, file_name = os.path.split(input_path)
    output_name = file_name.split(".", 1)[0] + ".asm"
    output_path = os.path.join(output_path, output_name)
    boot_code = get_boot_code(input_path)
    with open(output_path, "w") as output_file:
        output_file.write(boot_code)
    return output_path


def write_in_file(instructions: List[str], input_path: str) -> None:
    output_path = init_output_file(input_path)
    with open(output_path, "a") as output_file:
        for instruction in instructions:
            output_file.write(instruction + "\n")


def translate_file(vm_file_name: str) -> List[str]:
    instructions = parser.get_instructions(vm_file_name)
    asm_instructions = [
        translator.translate(instruction, vm_file_name) for instruction in instructions
    ]
    return asm_instructions


def translate(vm_file_or_directory_name: str) -> None:
    asm_instructions = []
    if os.path.isdir(vm_file_or_directory_name):
        for file_name in os.listdir(vm_file_or_directory_name):
            if file_name.endswith(".vm"):
                asm_instructions += translate_file(
                    vm_file_or_directory_name + "/" + file_name
                )
    else:
        asm_instructions = translate_file(vm_file_or_directory_name)
    write_in_file(asm_instructions, vm_file_or_directory_name)
