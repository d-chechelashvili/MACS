import os


def translate(instruction: str, vm_file: str) -> str:
    asm_translation = "// " + instruction + "\n"
    instruction_type = instruction.split()[0].upper()
    _, file_name = os.path.split(vm_file)
    asm_translation += FUNCTION_TABLE[instruction_type](instruction, file_name)
    return asm_translation


def translate_add(instruction: str, vm_file: str) -> str:
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "M=M+D" + "\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    return translation


def translate_sub(instruction: str, vm_file: str) -> str:
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "M=M-D" + "\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    return translation


def translate_neg(instruction: str, vm_file: str) -> str:
    translation = ("@SP" + "\n"
                   + "A=M-1" + "\n"
                   + "M=-M" + "\n")
    return translation


def translate_eq(instruction: str, vm_file: str) -> str:
    global logical_counter
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "D=M-D" + "\n"
        + "@TRUE" + str(logical_counter) + "\n"
        + "D;JEQ" + "\n"
        + "@SP" + "\n"
        + "A=M-1" + "\n"
        + "A=A-1" + "\n"
        + "M=0" + "\n"
        + "@ENDLOGICAL!" + str(logical_counter) + "\n"
        + "0;JMP" + "\n"
        + "(TRUE" + str(logical_counter) + ")\n"
        + "@SP" + "\n"
        + "A=M-1" + "\n"
        + "A=A-1" + "\n"
        + "M=-1" + "\n"
        + "(ENDLOGICAL!" + str(logical_counter) + ")\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    logical_counter += 1
    return translation


def translate_gt(instruction: str, vm_file: str) -> str:
    global logical_counter
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "D=M-D" + "\n"
        + "@TRUE" + str(logical_counter) + "\n"
        + "D;JGT" + "\n"
        + "@SP" + "\n"
        + "A=M-1" + "\n"
        + "A=A-1" + "\n"
        + "M=0" + "\n"
        + "@ENDLOGICAL!" + str(logical_counter) + "\n"
        + "0;JMP" + "\n"
        + "(TRUE" + str(logical_counter) + ")\n"
        + "@SP" + "\n"
        + "A=M-1" + "\n"
        + "A=A-1" + "\n"
        + "M=-1" + "\n"
        + "(ENDLOGICAL!" + str(logical_counter) + ")\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    logical_counter += 1
    return translation


def translate_lt(instruction: str, vm_file: str) -> str:
    global logical_counter
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "D=M-D" + "\n"
        + "@TRUE" + str(logical_counter) + "\n"
        + "D;JLT" + "\n"
        + "@SP" + "\n"
        + "A=M-1" + "\n"
        + "A=A-1" + "\n"
        + "M=0" + "\n"
        + "@ENDLOGICAL!" + str(logical_counter) + "\n"
        + "0;JMP" + "\n"
        + "(TRUE" + str(logical_counter) + ")\n"
        + "@SP" + "\n"
        + "A=M-1" + "\n"
        + "A=A-1" + "\n"
        + "M=-1" + "\n"
        + "(ENDLOGICAL!" + str(logical_counter) + ")\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    logical_counter += 1
    return translation


def translate_and(instruction: str, vm_file: str) -> str:
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "M=D&M" + "\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    return translation


def translate_or(instruction: str, vm_file: str) -> str:
    translation = (
        "@SP" + "\n"
        + "A=M-1" + "\n"
        + "D=M" + "\n"
        + "A=A-1" + "\n"
        + "M=D|M" + "\n"
        + "@SP" + "\n"
        + "M=M-1" + "\n"
    )
    return translation


def translate_not(instruction: str, vm_file: str) -> str:
    translation = ("@SP" + "\n"
                   + "A=M-1" + "\n"
                   + "M=!M" + "\n")
    return translation


def translate_push(instruction: str, vm_file: str) -> str:
    _, segment, number = instruction.split()
    translation = "UNDEFINED"
    if segment == "constant":
        translation = (
            "@" + number + "\n"
            + "D=A" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
        )
    elif segment in BASIC_SEGMENT_DICT:
        segment_addr = BASIC_SEGMENT_DICT[segment]
        translation = (
            "@" + number + "\n"
            + "D=A" + "\n"
            + "@" + segment_addr + "\n"
            + "A=M" + "\n"
            + "A=D+A" + "\n"
            + "D=M" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
        )
    elif segment == "temp":
        translation = (
            "@5" + "\n"
            + "D=A" + "\n"
            + "@" + number + "\n"
            + "A=D+A" + "\n"
            + "D=M" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
        )
    elif segment == "pointer":
        mapping = {"0": "THIS", "1": "THAT"}
        translation = (
            "@" + mapping[number] + "\n"
            + "D=M" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
        )
    elif segment == "static":
        translation = (
            "@" + vm_file + "." + number + "\n"
            + "D=M" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
        )

    return translation


def translate_pop(instruction: str, vm_file: str) -> str:
    _, segment, number = instruction.split()
    translation = "UNDEFINED"
    if segment in BASIC_SEGMENT_DICT:
        segment_addr = BASIC_SEGMENT_DICT[segment]
        translation = (
            "@" + number + "\n"
            + "D=A" + "\n"
            + "@" + segment_addr + "\n"
            + "A=M" + "\n"
            + "D=D+A" + "\n"
            + "@R14" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M-1" + "\n"
            + "A=M" + "\n"
            + "D=M" + "\n"
            + "@R14" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
        )
    elif segment == "temp":
        translation = (
            "@5" + "\n"
            + "D=A" + "\n"
            + "@" + number + "\n"
            + "D=D+A" + "\n"
            + "@R14" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M-1" + "\n"
            + "A=M" + "\n"
            + "D=M" + "\n"
            + "@R14" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
        )
    elif segment == "pointer":
        mapping = {"0": "THIS", "1": "THAT"}
        translation = (
            "@SP" + "\n"
            + "M=M-1" + "\n"
            + "A=M" + "\n"
            + "D=M" + "\n"
            + "@" + mapping[number] + "\n"
            + "M=D" + "\n"
        )
    elif segment == "static":
        translation = (
            "@SP" + "\n"
            + "M=M-1" + "\n"
            + "A=M" + "\n"
            + "D=M" + "\n"
            + "@" + vm_file + "." + number + "\n"
            + "M=D" + "\n"
        )
    return translation


logical_counter = 0

BASIC_SEGMENT_DICT = {"local": "LCL", "argument": "ARG", "this": "THIS", "that": "THAT"}

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
}
