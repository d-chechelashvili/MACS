from typing import List


def get_raw_instructions(file_name: str) -> List[str]:
    with open(file_name, "r") as file:
        file_content = file.read()
    raw_instructions = file_content.split("\n")
    return raw_instructions


def process_instructions(raw_instructions: List[str]) -> List[str]:
    instructions: List[str] = []
    for instruction in raw_instructions:
        instruction = instruction.strip()
        instruction = instruction.split("//", 1)[0]
        if instruction:
            instructions.append(instruction)
    return instructions


def get_instructions(file_name: str) -> List[str]:
    raw_instructions = get_raw_instructions(file_name)
    instructions = process_instructions(raw_instructions)
    return instructions
