class SymbolTable():
    def __init__(self):
        self.class_table = {}  # name - (type, kind, number)
        self.subroutine_table = {}  # name - (type, kind, number)
        self.next_index = {"static": 0, "this": 0, "argument": 0, "local": 0}

    def __contains__(self, item):
        return self.__search(item) is not None

    def start_subroutine(self):
        self.subroutine_table = {}
        self.next_index["argument"] = self.next_index["local"] = 0

    def define(self, name, type, kind):
        kind = kind if kind != "field" else "this"
        if kind in ["static", "this"]:
            self.__add_to_table(name, type, kind, self.class_table)
        elif kind in ["local", "argument"]:
            self.__add_to_table(name, type, kind, self.subroutine_table)
        else:
            print("HMM?:", name, type, kind)

    def __add_to_table(self, name, type, kind, table):
        table[name] = (type, kind, self.var_count(kind))
        self.next_index[kind] += 1

    def var_count(self, kind):
        return self.next_index[kind]

    def type_of(self, name):
        return self.__get_argument(name, 0)

    def kind_of(self, name):
        return self.__get_argument(name, 1)

    def index_of(self, name):
        return self.__get_argument(name, 2)

    def print_entry(self, name):
        print("TOKEN:", name,
              "TYPE:", self.type_of(name),
              "KIND:", self.kind_of(name),
              "INDEX:", self.index_of(name))

    def __get_argument(self, name, index):
        entry = self.__search(name)
        if entry: return entry[index]
        return None

    def __search(self, name):
        if name in self.subroutine_table:
            return self.subroutine_table[name]
        elif name in self.class_table:
            return self.class_table[name]
        return None
