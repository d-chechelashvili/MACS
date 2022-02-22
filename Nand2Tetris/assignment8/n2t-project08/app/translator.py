import os

from app.translator_I import (
    translate_add,
    translate_and,
    translate_eq,
    translate_gt,
    translate_lt,
    translate_neg,
    translate_not,
    translate_or,
    translate_pop,
    translate_push,
    translate_sub,
)
from app.translator_II import (
    translate_call,
    translate_conditional_branch,
    translate_function,
    translate_label,
    translate_return,
    translate_unconditional_branch,
)


def translate(instruction: str, vm_file: str) -> str:
    asm_translation = "// " + instruction + "\n"
    instruction_type = instruction.split()[0].upper()
    _, file_name = os.path.split(vm_file)
    asm_translation += FUNCTION_TABLE[instruction_type](instruction, file_name)
    return asm_translation


FUNCTION_TABLE = {
    "ADD": translate_add,
    "SUB": translate_sub,
    "NEG": translate_neg,
    "EQ": translate_eq,
    "GT": translate_gt,
    "LT": translate_lt,
    "AND": translate_and,
    "OR": translate_or,
    "NOT": translate_not,
    "PUSH": translate_push,
    "POP": translate_pop,
    "LABEL": translate_label,
    "GOTO": translate_unconditional_branch,
    "IF-GOTO": translate_conditional_branch,
    "FUNCTION": translate_function,
    "CALL": translate_call,
    "RETURN": translate_return,
}
