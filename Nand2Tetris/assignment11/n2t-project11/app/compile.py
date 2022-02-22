# Your code starts here:
import os

from app import utils
from app.compilation_engine import CompilationEngine
from app.jack_tokenizer import JackTokenizer
from app.symbol_table import SymbolTable


def compile_file(jack_file_name: str) -> None:
    compilation_engine = CompilationEngine(jack_file_name)
    compilation_engine.compile_class()


def compile(jack_file_or_directory_name: str) -> None:
    if os.path.isdir(jack_file_or_directory_name):
        for file_name in os.listdir(jack_file_or_directory_name):
            if file_name.endswith(".jack"):
                compile_file(jack_file_or_directory_name + "/" + file_name)
    else:
        compile_file(jack_file_or_directory_name)
