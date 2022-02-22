from app.constants import comp_table, dest_table, jump_table
from app.SymbolTable import SymbolTable


def translate_A(instruction: str, symbol_table: SymbolTable) -> str:
    address = instruction[1:]
    binary_representation = "0"
    if not address.isdigit():
        address = symbol_table.get_address(address)
    binary_address = bin(int(address)).replace("0b", "").zfill(15)
    binary_representation += binary_address
    return binary_representation


def translate_dest(dest: str) -> str:
    return dest_table.get(dest, "")


def translate_comp(comp: str) -> str:
    binary_comp = "1" if "M" in comp else "0"
    binary_comp += comp_table.get(comp, "")
    return binary_comp


def translate_jump(jump: str) -> str:
    return jump_table.get(jump, "")


def translate_c(instruction: str) -> str:
    binary_representation = "111"
    dest = comp = jump = ""

    if "=" not in instruction:
        comp, jump = instruction.split(";", 1)
    else:
        dest, comp_jump = instruction.split("=", 1)
        if ";" not in comp_jump:
            comp = comp_jump
        else:
            comp, jump = instruction.split(";", 1)

    binary_representation += (
        translate_comp(comp) + translate_dest(dest) + translate_jump(jump)
    )
    return binary_representation


def translate(instruction: str, symbol_table: SymbolTable) -> str:
    binary_representation = (
        translate_A(instruction, symbol_table)
        if instruction[0] == "@"
        else translate_c(instruction)
    )
    return binary_representation
