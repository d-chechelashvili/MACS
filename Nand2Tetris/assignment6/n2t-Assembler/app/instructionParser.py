import re
from typing import List

from app.SymbolTable import SymbolTable


def get_raw_instructions(file_name: str) -> List[str]:
    with open(file_name, "r") as file:
        file_content = file.read()
    raw_instructions = file_content.split("\n")
    return raw_instructions


def process_instructions(
    raw_instructions: List[str], symbol_table: SymbolTable
) -> List[str]:
    instructions: List[str] = []
    labels: List[str] = []
    find_label_address = False
    for instruction in raw_instructions:
        instruction = "".join(instruction.split())  # remove spaces
        instruction = instruction.split("//", 1)[0]  # remove comments
        # if there is a label on this line
        # when we get next real instruction we add this label to symbolTable
        label = re.search(r"\((.*?)\)", instruction)
        if label:
            labels.append(label.group(1))
            instruction = ""
            find_label_address = True
        if instruction:
            if find_label_address:
                for label in labels:
                    symbol_table.add_entry(label, len(instructions))
                labels = []
                find_label_address = False
            instructions.append(instruction)
    return instructions


def get_instructions(file_name: str, symbol_table: SymbolTable) -> List[str]:
    raw_instructions = get_raw_instructions(file_name)
    instructions = process_instructions(raw_instructions, symbol_table)
    return instructions
