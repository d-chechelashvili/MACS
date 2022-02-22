logical_counter = 0

BASIC_SEGMENT_MAPPING = {
    "local": "LCL",
    "argument": "ARG",
    "this": "THIS",
    "that": "THAT",
}

POINTER_MAPPING = {"0": "THIS", "1": "THAT"}


def translate_add(instruction: str, vm_file: str) -> str:
    return translate_operation("+")


def translate_sub(instruction: str, vm_file: str) -> str:
    return translate_operation("-")


def translate_neg(instruction: str, vm_file: str) -> str:
    return translate_negation("-")


def translate_eq(instruction: str, vm_file: str) -> str:
    return translate_comparison("EQ")


def translate_gt(instruction: str, vm_file: str) -> str:
    return translate_comparison("GT")


def translate_lt(instruction: str, vm_file: str) -> str:
    return translate_comparison("LT")


def translate_and(instruction: str, vm_file: str) -> str:
    return translate_operation("&")


def translate_or(instruction: str, vm_file: str) -> str:
    return translate_operation("|")


def translate_not(instruction: str, vm_file: str) -> str:
    return translate_negation("!")


def translate_operation(operator: str) -> str:
    translation = (
            "@SP" + "\n"
            + "A=M-1" + "\n"
            + "D=M" + "\n"
            + "A=A-1" + "\n"
            + "M=M" + operator + "D" + "\n"
            + "@SP" + "\n"
            + "M=M-1" + "\n"
    )
    return translation


def translate_comparison(operator: str) -> str:
    global logical_counter
    translation = (
            "@SP" + "\n"
            + "A=M-1" + "\n"
            + "D=M" + "\n"
            + "A=A-1" + "\n"
            + "D=M-D" + "\n"
            + "@TRUE" + str(logical_counter) + "\n"
            + "D;J" + operator + "\n"
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


def translate_negation(operator: str) -> str:
    translation = (
            "@SP" + "\n"
            + "A=M-1" + "\n"
            + "M=" + operator + "M" + "\n"
    )
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
    elif segment in BASIC_SEGMENT_MAPPING:
        segment_addr = BASIC_SEGMENT_MAPPING[segment]
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
        translation = (
                "@" + POINTER_MAPPING[number] + "\n"
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
    if segment in BASIC_SEGMENT_MAPPING:
        segment_addr = BASIC_SEGMENT_MAPPING[segment]
        translation = (
                "@" + number + "\n"
                + "D=A" + "\n"
                + "@" + segment_addr + "\n"
                + "A=M" + "\n"
                + "D=D+A" + "\n"
                + "@R13" + "\n"
                + "M=D" + "\n"
                + "@SP" + "\n"
                + "M=M-1" + "\n"
                + "A=M" + "\n"
                + "D=M" + "\n"
                + "@R13" + "\n"
                + "A=M" + "\n"
                + "M=D" + "\n"
        )
    elif segment == "temp":
        translation = (
                "@5" + "\n"
                + "D=A" + "\n"
                + "@" + number + "\n"
                + "D=D+A" + "\n"
                + "@R13" + "\n"
                + "M=D" + "\n"
                + "@SP" + "\n"
                + "M=M-1" + "\n"
                + "A=M" + "\n"
                + "D=M" + "\n"
                + "@R13" + "\n"
                + "A=M" + "\n"
                + "M=D" + "\n"
        )
    elif segment == "pointer":
        translation = (
                "@SP" + "\n"
                + "M=M-1" + "\n"
                + "A=M" + "\n"
                + "D=M" + "\n"
                + "@" + POINTER_MAPPING[number] + "\n"
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
