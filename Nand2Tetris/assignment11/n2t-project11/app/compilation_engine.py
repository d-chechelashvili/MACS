from app.VM_writer import VMWriter
from app.jack_tokenizer import JackTokenizer
from app.symbol_table import SymbolTable


class CompilationEngine:

    def __init__(self, input_path: str):
        self.__tokenizer = JackTokenizer(input_path)
        self.__vm_writer = VMWriter(input_path)
        self.__symbol_table = SymbolTable()
        self.__counter = 0
        self.class_name = ""
        self.__statement_functions = {
            "do": self.compile_do,
            "let": self.compile_let,
            "while": self.compile_while,
            "return": self.compile_return,
            "if": self.compile_if
        }
        self.__mappings = {
            "+": "add",
            "-": "sub",
            "*": "Math.multiply",
            "/": "Math.divide",
            "&amp;": "and",
            "|": "or",
            "&lt;": "lt",
            "&gt;": "gt",
            "=": "eq"
        }

    # contract - after each method pointer is standing just after that method

    def advance(self):
        self.__tokenizer.advance()

    def type(self):
        return self.__tokenizer.token_type()

    def token(self):
        return self.__tokenizer.token()

    def define(self, name, type, kind):
        self.__symbol_table.define(name, type, kind)

    # ---------------------------------------------------------------

    def compile_class(self):
        self.advance()  # class
        self.advance()  # className
        self.class_name = self.token()
        self.advance()  # {
        self.advance()  # classVarDec/subroutineDec/}
        while self.token() in ["static", "field"]:
            self.compile_class_var_dec()

        while self.token() in ["constructor", "function", "method"]:
            self.compile_subroutine_dec()
        # I am standing on } so I dont really have to do anything

    def compile_class_var_dec(self):
        # I am standing on variable declaration type (static/field)
        kind = self.token()
        self.advance()  # type
        type = self.token()
        while True:
            self.advance()  # varName
            self.define(self.token(), type, kind)
            self.advance()  # ,/;
            if self.token() != ",":
                break
        self.advance()

    def compile_subroutine_dec(self):
        # I am standing on function/method/constructor
        self.__symbol_table.start_subroutine()
        self.subroutine_type = self.token()
        if self.subroutine_type == "method":
            self.define("this", self.class_name, "argument")

        self.advance()  # type
        self.advance()  # subroutineName
        self.subroutine_name = self.token()
        self.advance()  # (
        self.advance()  # parameterList / )
        self.compile_parameter_list()

        # I am standing on )
        self.advance()  # subroutineBody
        self.compile_subroutine_body()

    def compile_subroutine_body(self):
        # I am standing on {
        self.advance()  # varDec/statements
        local_var_num = 0
        while self.token() == "var":
            local_var_num += self.compile_var_dec()

        self.__init_subroutine(local_var_num)
        # I am standing on statements
        self.compile_statements()
        # I am standing on }
        self.advance()

    def __init_subroutine(self, local_var_num):
        self.__vm_writer.write_function(self.class_name + "." + self.subroutine_name, local_var_num)
        if self.subroutine_type == "method":
            self.__vm_writer.write_push("argument", 0)
            self.__vm_writer.write_pop("pointer", 0)
        elif self.subroutine_type == "constructor":
            self.__vm_writer.write_push("constant", self.__symbol_table.var_count('this'))
            self.__vm_writer.write_call("Memory.alloc", 1)
            self.__vm_writer.write_pop("pointer", 0)

    def compile_parameter_list(self):
        while self.token() != ")":
            # I am standing on type of parameter
            type = self.type()
            self.advance()  # varName
            self.define(self.token(), type, "argument")
            self.advance()  # ,/)
            if self.token() == ",":
                self.advance()

    def compile_var_dec(self):
        var_num = 0
        # I am standing on var
        self.advance()  # type
        type = self.token()
        while True:
            var_num += 1
            self.advance()  # varName
            var_name = self.token()
            self.__symbol_table.define(var_name, type, 'local')
            self.advance()  # ,/;
            if self.token() == ";":
                break

        # I am standing on ; so I should advance once more according to contract
        self.advance()
        return var_num

    def compile_statements(self):
        # I am standing on statement do/let/while/return/if
        while self.token() in ["do", "let", "while", "return", "if"]:
            self.__statement_functions[self.token()]()

    def compile_do(self):
        # I am standing on do
        self.advance()  # subroutineCall
        self.__make_subroutine_call()
        self.__vm_writer.write_pop('temp', 0)  # pop return value
        # I am standing on ;
        self.advance()

    def __make_subroutine_call(self):
        # I am standing on subroutineName/className/varName
        name = self.token()
        self.advance()  # (/.
        args_num = 1
        class_name = self.class_name
        if self.token() == ".":
            if name not in self.__symbol_table:
                args_num = 0
                class_name = name
            else:
                class_name = self.__symbol_table.type_of(name)
                self.__vm_writer.write_push(self.__symbol_table.kind_of(name),
                                            self.__symbol_table.index_of(name))

            self.advance()  # subroutineName
            name = self.token()
            self.advance()  # (
        else:
            self.__vm_writer.write_push('pointer', 0)

        subroutine_name = class_name + "." + name
        self.advance()  # expressionList
        args_num += self.compile_expression_list()
        self.__vm_writer.write_call(subroutine_name, args_num)
        # I am standing on )
        self.advance()

    def compile_let(self):
        # I am standing on let
        self.advance()  # varName
        var_name = self.token()
        self.advance()  # [/=
        if self.token() == "[":
            self.advance()  # expression
            self.__vm_writer.write_push(self.__symbol_table.kind_of(var_name),
                                        self.__symbol_table.index_of(var_name))
            self.compile_expression()
            self.__vm_writer.write_arithmetic("add")
            # I am standing on ]
            self.advance()  # =
            self.advance()  # expression
            self.compile_expression()
            # array access code
            self.__vm_writer.write_pop("temp", 0)
            self.__vm_writer.write_pop("pointer", 1)
            self.__vm_writer.write_push("temp", 0)
            self.__vm_writer.write_pop("that", 0)
            # I am standing on ;
        else:
            self.advance()  # expression
            self.compile_expression()
            self.__vm_writer.write_pop(self.__symbol_table.kind_of(var_name),
                                       self.__symbol_table.index_of(var_name))
            # I am standing on ;

        self.advance()

    def compile_while(self):
        L1 = "L1" + str(self.__counter)
        L2 = "L2" + str(self.__counter)
        self.__counter += 1
        # I am standing on while
        self.advance()  # (
        self.advance()  # expression
        self.__vm_writer.write_label(L1)
        self.compile_expression()
        self.__vm_writer.write_arithmetic('not')
        self.__vm_writer.write_if(L2)
        # I am standing on )
        self.advance()  # {
        self.advance()  # statements
        self.compile_statements()
        self.__vm_writer.write_go_to(L1)
        self.__vm_writer.write_label(L2)
        # I am standing on }
        self.advance()  # after while

    def compile_return(self):
        # I am standing on return
        self.advance()  # expression/;
        if self.token() == ";":
            self.__vm_writer.write_push("constant", 0)  # push 0 when returning nothing
        else:
            self.compile_expression()
        self.__vm_writer.write_return()
        # I am standing on ;
        self.advance()  # advance after return

    def compile_if(self):
        L1 = "L1" + str(self.__counter)
        L2 = "L2" + str(self.__counter)
        self.__counter += 1
        # I am standing on if
        self.advance()  # (
        self.advance()  # expression
        self.compile_expression()
        self.__vm_writer.write_arithmetic('not')
        # I am standing on )
        self.advance()  # {
        self.advance()  # statements
        self.__vm_writer.write_if(L1)
        self.compile_statements()
        self.__vm_writer.write_go_to(L2)
        # I am standing on }
        self.advance()  # else/ end
        self.__vm_writer.write_label(L1)
        if self.token() == "else":
            self.advance()  # {
            self.advance()  # statements
            self.compile_statements()
            # I am standing on }
            self.advance()  # after if
        self.__vm_writer.write_label(L2)

    def compile_expression_list(self):
        args_num = 0
        while self.token() != ")":
            args_num += 1
            # I am standing on expression
            self.compile_expression()
            # I am standing on , or )
            if self.token() == ",":
                self.advance()  # expression
        return args_num

    def compile_expression(self):
        # I am standing on term
        self.compile_term()
        # I am standing on either operator or after term
        while self.token() in ["+", "-", "*", "/", "&amp;", "|", "&lt;", "&gt;", "="]:
            operator = self.token()
            self.advance()  # term
            self.compile_term()
            command = self.__mappings[operator]
            if operator in ["*", "/"]:
                self.__vm_writer.write_call(command, 2)  # translate to Math.multiply/divide
            else:
                self.__vm_writer.write_arithmetic(command)

    def compile_term(self):
        # I am standing on either integerConstant or stringConstant or keywordConstant
        # or varname or subroutineCall or ( or unaryOp
        token = self.token()
        if self.token() == "(":
            self.advance()  # expression
            self.compile_expression()
            # I am standing on )
            self.advance()  # advance after term
        elif self.token() == "-" or self.token() == "~":
            self.advance()  # term
            self.compile_term()
            arithmetic = "neg" if token == "-" else "not"
            self.__vm_writer.write_arithmetic(arithmetic)
        elif self.type() != "identifier":
            self.__compile_term_constants()
        else:
            self.__compile_term_identifier()

    def __compile_term_constants(self):
        # I am standing on integerConstant/stringConstant/keywordConstant
        if self.type() == "integerConstant":
            self.__vm_writer.write_push("constant", self.token())
        elif self.type() == "keyword":
            # its either this or true or false or null
            if self.token() == "this":
                self.__vm_writer.write_push("pointer", 0)
            elif self.token() == "true":
                self.__vm_writer.write_push("constant", 0)
                self.__vm_writer.write_arithmetic("not")
            else:
                self.__vm_writer.write_push("constant", 0)
        else:
            # its string constant
            # create length argument to pass to string.new
            length = len(self.token())
            self.__vm_writer.write_push("constant", length)
            # create empty string
            self.__vm_writer.write_call("String.new", 1)
            for char in self.token():
                # push char integer value on stack
                self.__vm_writer.write_push("constant", ord(char))
                # invoke append function with char and pointer to string object arguments
                self.__vm_writer.write_call("String.appendChar", 2)

        self.advance()  # advance after term

    def __compile_term_identifier(self):
        # I am standing on either varName, or varName[expresion] or subroutineCall
        self.advance()  # after/[/(/.
        if self.token() == "[":
            var_name = self.__tokenizer.prev()[0]
            self.__vm_writer.write_push(self.__symbol_table.kind_of(var_name),
                                        self.__symbol_table.index_of(var_name))
            self.advance()  # expression
            self.compile_expression()
            self.__vm_writer.write_arithmetic("add")
            self.__vm_writer.write_pop("pointer", 1)
            self.__vm_writer.write_push("that", 0)
            # I am standing on ]
            self.advance()
        elif self.token() == "(" or self.token() == ".":
            self.__tokenizer.rewind()  # go back on subroutineName
            self.__make_subroutine_call()
        else:
            # I was standing on varName
            var_name = self.__tokenizer.prev()[0]
            segment = self.__symbol_table.kind_of(var_name)
            index = self.__symbol_table.index_of(var_name)
            self.__vm_writer.write_push(segment, index)
