from app import utils
from app.jack_tokenizer import JackTokenizer


class CompilationEngine:

    def __init__(self, input_path: str):
        self.__output_path = utils.init_output_file(input_path, tokenizer_log=False)
        self.__tokenizer = JackTokenizer(input_path, self.__output_path)
        self.__statement_functions = {
            'do': self.compile_do,
            'let': self.compile_let,
            'while': self.compile_while,
            'return': self.compile_return,
            'if': self.compile_if,
        }

    def compile_class(self):
        utils.write_in_file("<class>\n", self.__output_path)
        self.__tokenizer.advance()  # class
        self.__tokenizer.advance()  # className
        self.__tokenizer.advance()  # {
        self.__tokenizer.advance(log=False)  # classvardec/subroutinedec/}
        while self.__tokenizer.token() in ["static", "field"]:
            self.compile_class_var_dec()
            self.__tokenizer.advance(log=False)  # classvardec/subroutinedec/}

        while self.__tokenizer.token() in ["constructor", "function", "method"]:
            self.compile_subroutine_dec()
            self.__tokenizer.advance(log=False)  # subroutinedec/}

        # assert self.__tokenizer.token() == "}"
        self.__tokenizer.log_token()
        utils.write_in_file("</class>\n", self.__output_path)

    def compile_class_var_dec(self):
        utils.write_in_file("<classVarDec>\n", self.__output_path)
        self.__tokenizer.log_token()  # static/field
        self.__tokenizer.advance()  # type
        self.__tokenizer.advance()  # varName

        self.__tokenizer.advance()  # ,/;
        while self.__tokenizer.token() == ",":
            self.__tokenizer.advance()  # varName
            self.__tokenizer.advance()  # ,/;
        utils.write_in_file("</classVarDec>\n", self.__output_path)

    def compile_subroutine_dec(self):
        utils.write_in_file("<subroutineDec>\n", self.__output_path)
        self.__tokenizer.log_token()  # constructor/func/method
        self.__tokenizer.advance()  # void/type
        self.__tokenizer.advance()  # subroutineName
        self.__tokenizer.advance()  # (

        self.__tokenizer.advance(log=False)  # parameterList / )
        self.compile_parameter_list()

        self.__tokenizer.log_token()  # )

        self.__tokenizer.advance(log=False)  # subroutineBody
        self.compile_subroutine_body()
        utils.write_in_file("</subroutineDec>\n", self.__output_path)

    def compile_subroutine_body(self):
        utils.write_in_file("<subroutineBody>\n", self.__output_path)
        self.__tokenizer.log_token()  # {

        self.__tokenizer.advance(log=False)  # varDec/statements
        while self.__tokenizer.token() == "var":
            self.compile_var_dec()
            self.__tokenizer.advance(log=False)  # varDec/statements

        self.compile_statements()

        self.__tokenizer.log_token()  # }
        utils.write_in_file("</subroutineBody>\n", self.__output_path)

    def compile_parameter_list(self):
        utils.write_in_file("<parameterList>\n", self.__output_path)

        while self.__tokenizer.token() != ")":
            self.__tokenizer.log_token()  # type
            self.__tokenizer.advance()  # varName
            self.__tokenizer.advance(log=False)  # ,/)
            if self.__tokenizer.token() == ",":
                self.__tokenizer.log_token()
                self.__tokenizer.advance(log=False)

        utils.write_in_file("</parameterList>\n", self.__output_path)

    def compile_var_dec(self):
        utils.write_in_file("<varDec>\n", self.__output_path)

        self.__tokenizer.log_token()  # var
        self.__tokenizer.advance()  # type
        self.__tokenizer.advance()  # varName

        self.__tokenizer.advance()  # ,/;
        while self.__tokenizer.token() == ",":
            self.__tokenizer.advance()  # varName
            self.__tokenizer.advance()  # ,/;

        utils.write_in_file("</varDec>\n", self.__output_path)

    def compile_statements(self):
        utils.write_in_file("<statements>\n", self.__output_path)

        while self.__tokenizer.token() in ["do", "let", "while", "return", "if"]:
            self.__statement_functions[self.__tokenizer.token()]()

        utils.write_in_file("</statements>\n", self.__output_path)

    def compile_do(self):
        utils.write_in_file("<doStatement>\n", self.__output_path)
        self.__tokenizer.log_token()  # do
        self.__make_subroutine_call()
        self.__tokenizer.log_token()  # ;
        self.__tokenizer.advance(log=False)
        utils.write_in_file("</doStatement>\n", self.__output_path)

    def __make_subroutine_call(self, log_name=True):
        if log_name:
            self.__tokenizer.advance()  # subroutineName /className/varName
            self.__tokenizer.advance()  # ( / .
        if self.__tokenizer.token() == ".":
            self.__tokenizer.advance()  # subroutineName
            self.__tokenizer.advance()  # (
        self.__tokenizer.advance(log=False)  # expressionList/)
        self.compile_expression_list()
        self.__tokenizer.log_token()  # )
        if log_name: self.__tokenizer.advance(log=False)

    def compile_let(self):
        utils.write_in_file("<letStatement>\n", self.__output_path)
        self.__tokenizer.log_token()  # let
        self.__tokenizer.advance()  # varName
        self.__tokenizer.advance()  # [/=
        if self.__tokenizer.token() == "[":
            self.__tokenizer.advance(log=False)
            self.compile_expression()
            self.__tokenizer.log_token()  # ]
            self.__tokenizer.advance()  # =
        self.__tokenizer.advance(log=False)  # expression
        self.compile_expression()
        self.__tokenizer.log_token()  # ;
        self.__tokenizer.advance(log=False)
        utils.write_in_file("</letStatement>\n", self.__output_path)

    def compile_while(self):
        utils.write_in_file("<whileStatement>\n", self.__output_path)
        self.__tokenizer.log_token()  # while
        self.__tokenizer.advance()  # (
        self.__tokenizer.advance(log=False)
        self.compile_expression()
        self.__tokenizer.log_token()  # )
        self.__tokenizer.advance()  # {
        self.__tokenizer.advance(log=False)
        self.compile_statements()
        self.__tokenizer.log_token()  # }
        self.__tokenizer.advance(log=False)
        utils.write_in_file("</whileStatement>\n", self.__output_path)

    def compile_return(self):
        utils.write_in_file("<returnStatement>\n", self.__output_path)
        self.__tokenizer.log_token()  # return
        self.__tokenizer.advance(log=False)  # expression/;
        if self.__tokenizer.token() != ";":
            self.compile_expression()
        self.__tokenizer.log_token()  # ;
        self.__tokenizer.advance(log=False)
        utils.write_in_file("</returnStatement>\n", self.__output_path)

    def compile_if(self):
        utils.write_in_file("<ifStatement>\n", self.__output_path)
        self.__tokenizer.log_token()  # if
        self.__tokenizer.advance()  # (
        self.__tokenizer.advance(log=False)
        self.compile_expression()
        self.__tokenizer.log_token()  # )
        self.__tokenizer.advance()  # {
        self.__tokenizer.advance(log=False)
        self.compile_statements()
        self.__tokenizer.log_token()  # }
        self.__tokenizer.advance(log=False)  # else?
        if self.__tokenizer.token() == "else":
            self.__tokenizer.log_token()
            self.__tokenizer.advance()  # {
            self.__tokenizer.advance(log=False)
            self.compile_statements()
            self.__tokenizer.log_token()  # }
            self.__tokenizer.advance(log=False)
        utils.write_in_file("</ifStatement>\n", self.__output_path)

    def compile_expression(self):
        utils.write_in_file("<expression>\n", self.__output_path)
        self.compile_term()
        while self.__tokenizer.token() in ["+", "-", "*", "/", "&amp;", "|", "&lt;", "&gt;", "="]:
            self.__tokenizer.log_token()  # op
            self.__tokenizer.advance(log=False)  # term
            self.compile_term()
        utils.write_in_file("</expression>\n", self.__output_path)

    def compile_term(self):
        utils.write_in_file("<term>\n", self.__output_path)

        if self.__tokenizer.token() == "(":
            self.__tokenizer.log_token()  # (
            self.__tokenizer.advance(log=False)
            self.compile_expression()
            self.__tokenizer.log_token()  # )
            self.__tokenizer.advance(log=False)  # end
        elif self.__tokenizer.token() == "-" or self.__tokenizer.token() == "~":
            self.__tokenizer.log_token()  # -/~
            self.__tokenizer.advance(log=False)  # term
            self.compile_term()
        elif self.__tokenizer.token_type() != "identifier":
            self.__tokenizer.log_token()  # integerConstant/stringConstant/keywordConstant
            self.__tokenizer.advance(log=False)  # end
        else:
            # we have either varName, or varName[expresion] or subroutineCall
            self.__tokenizer.log_token()  # varName, subroutineName, className
            self.__tokenizer.advance(log=False)  # end / [ / (

            if self.__tokenizer.token() == "[":
                self.__tokenizer.log_token()  # [
                self.__tokenizer.advance(log=False)  # expression
                self.compile_expression()
                self.__tokenizer.log_token()  # ]
                self.__tokenizer.advance(log=False)  # end
            elif self.__tokenizer.token() == "(" or self.__tokenizer.token() == ".":
                self.__tokenizer.log_token()
                self.__make_subroutine_call(log_name=False)
                self.__tokenizer.advance(log=False)  # end

        utils.write_in_file("</term>\n", self.__output_path)

    def compile_expression_list(self):
        utils.write_in_file("<expressionList>\n", self.__output_path)

        while self.__tokenizer.token() != ")":
            self.compile_expression()
            if self.__tokenizer.token() == ",":
                self.__tokenizer.log_token()  # ,
                self.__tokenizer.advance(log=False)  # expression

        utils.write_in_file("</expressionList>\n", self.__output_path)
