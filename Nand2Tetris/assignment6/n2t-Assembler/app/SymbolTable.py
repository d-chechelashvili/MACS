from typing import Dict


class SymbolTable:
    def __init__(self) -> None:
        self.__symbol_table: Dict[str, int] = {}
        self.__last_var_address: int = 15
        self.__add_predefined_symbols()

    def add_entry(self, symbol: str, address: int) -> None:
        self.__symbol_table[symbol] = address

    def __contains__(self, item: str) -> bool:
        return item in self.__symbol_table

    def get_address(self, symbol: str) -> int:
        if symbol not in self.__symbol_table:
            self.__last_var_address += 1
            self.add_entry(symbol, self.__last_var_address)
        return self.__symbol_table[symbol]

    def __add_predefined_symbols(self) -> None:
        for i in range(16):
            self.add_entry("R" + str(i), i)
        self.add_entry("SCREEN", 16384)
        self.add_entry("KBD", 24576)
        self.add_entry("SP", 0)
        self.add_entry("LCL", 1)
        self.add_entry("ARG", 2)
        self.add_entry("THIS", 3)
        self.add_entry("THAT", 4)
