import os

from app.jack_tokenizer import JackTokenizer
from app.compilation_engine import CompilationEngine


def analyze_file(jack_file_name: str) -> None:
    compilation_engine = CompilationEngine(jack_file_name)
    compilation_engine.compile_class()


def analyze(jack_file_or_directory_name: str) -> None:
    if os.path.isdir(jack_file_or_directory_name):
        for file_name in os.listdir(jack_file_or_directory_name):
            if file_name.endswith(".jack"):
                analyze_file(jack_file_or_directory_name + "/" + file_name)
    else:
        analyze_file(jack_file_or_directory_name)
