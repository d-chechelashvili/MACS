import os.path
from collections import defaultdict
from typing import Dict

from app.translator_I import translate_pop, translate_push

function_calls: Dict[str, int] = defaultdict(int)


def get_boot_code(input_path: str) -> str:
    boot_code = "//SP = 256\n"
    if os.path.isdir(input_path):
        boot_code += (
                "@256" + "\n"
                + "D=A" + "\n"
                + "@SP" + "\n"
                + "M=D" + "\n"
                + "\n//call Sys.init 0\n"
                + translate_call("call Sys.init 0", "Sys")
        )
    return boot_code


def translate_label(instruction: str, vm_file: str) -> str:
    _, label = instruction.split()
    translation = "(" + label + ")\n"
    return translation


def translate_unconditional_branch(instruction: str, vm_file: str) -> str:
    _, label = instruction.split()
    translation = (
            "@" + label + "\n"
            + "0;JMP" + "\n"
    )
    return translation


def translate_conditional_branch(instruction: str, vm_file: str) -> str:
    _, label = instruction.split()
    translation = (
            "@SP" + "\n"
            + "M=M-1" + "\n"
            + "A=M" + "\n"
            + "D=M" + "\n"
            + "@" + label + "\n"
            + "D;JNE" + "\n"
    )
    return translation


def translate_function(instruction: str, vm_file: str) -> str:
    _, function_name, variable_number = instruction.split()
    translation = translate_label("label " + function_name, vm_file)
    for i in range(int(variable_number)):
        translation += translate_push("push constant 0", vm_file)
    return translation


def translate_call(instruction: str, vm_file: str) -> str:
    _, function_name, argument_number = instruction.split()
    returnLabel = function_name + "$ret." + str(function_calls[function_name])
    translation = (
            "@" + returnLabel + "\n"
            + "D=A" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
            + __push("LCL")
            + __push("ARG")
            + __push("THIS")
            + __push("THAT")
            + "@" + argument_number + "\n"
            + "D=A" + "\n"
            + "@5" + "\n"
            + "D=D+A" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "D=A-D" + "\n"
            + "@ARG" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "D=M" + "\n"
            + "@LCL" + "\n"
            + "M=D" + "\n"
            + translate_unconditional_branch("goto " + function_name, vm_file)
            + translate_label("label " + returnLabel, vm_file)
    )
    function_calls[function_name] += 1
    return translation


def translate_return(instruction: str, vm_file: str) -> str:
    translation = (
            __load("R14", "5")
            + translate_pop("pop argument 0", vm_file)
            + "@ARG" + "\n"
            + "D=M+1" + "\n"
            + "@SP" + "\n"
            + "M=D" + "\n"
            + __load("THAT", "1")
            + __load("THIS", "2")
            + __load("ARG", "3")
            + __load("LCL", "4")
            + "@R14" + "\n"
            + "A=M" + "\n"
            + "0;JMP" + "\n"
    )
    return translation


# translates push address in stack. e.g. "push LCL"
def __push(address: str) -> str:
    translation = (
            "@" + address + "\n"
            + "D=M" + "\n"
            + "@SP" + "\n"
            + "A=M" + "\n"
            + "M=D" + "\n"
            + "@SP" + "\n"
            + "M=M+1" + "\n"
    )
    return translation


# translates destination = *(endFrame-offset). e.g. "THAT=*(endFrame-1)"
def __load(destination: str, offset: str) -> str:
    translation = (
            "@LCL" + "\n"
            + "D=M" + "\n"
            + "@" + offset + "\n"
            + "A=D-A" + "\n"
            + "D=M" + "\n"
            + "@" + destination + "\n"
            + "M=D" + "\n"
    )
    return translation
