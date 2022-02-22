import re
from typing import Tuple

from app import utils
from rop import Fail, Maybe, Railway, Some, Split


def parse_keyword(string: str) -> Maybe[Tuple[str, str]]:
    matchObj = re.match(
        r"class |constructor |function |method |field |static |var |int|char|boolean|void|true|false|null|this|let |do |if|else|while|return",
        string,
    )
    if matchObj:
        return Some((matchObj.group(), "keyword"))
    return Fail()


def parse_symbol(string: str) -> Maybe[Tuple[str, str]]:
    matchObj = re.match(r"[{}()[\].,;+\-*/&|<>=~]", string)
    if matchObj:
        return Some((matchObj.group(), "symbol"))
    return Fail()


def parse_integer_constant(string: str) -> Maybe[Tuple[str, str]]:
    matchObj = re.match(r"\d+", string)
    if matchObj:
        return Some((matchObj.group(), "integerConstant"))
    return Fail()


def parse_string_constant(string: str) -> Maybe[Tuple[str, str]]:
    matchObj = re.match(r"\"(.*?)\"", string)
    if matchObj:
        return Some((matchObj.group(), "stringConstant"))
    return Fail()


def parse_identifier(string: str) -> Maybe[Tuple[str, str]]:
    matchObj = re.match(r"(_|[a-z]|[A-Z])\w*", string)
    if matchObj:
        return Some((matchObj.group(), "identifier"))
    return Fail()


class JackTokenizer:
    def __init__(self, input_path: str, output_path: str) -> None:
        self.__input_path = input_path
        self.__output_path = output_path
        self.__code_string = self.__process_file()
        self.__token = ""
        self.__token_type = None

    def has_more_tokens(self) -> bool:
        return len(self.__code_string) > 0

    @property
    def __parse_next_token(self):
        return Railway(
            Split(
                parse_integer_constant,
                parse_keyword,
                parse_identifier,
                parse_symbol,
                parse_string_constant,
            ),
        )

    def advance(self, log=True) -> None:
        assert self.has_more_tokens()
        token_info = self.__parse_next_token(self.__code_string)
        assert isinstance(token_info, Some)

        self.__token, self.__token_type = token_info.value
        assert self.__code_string.startswith(self.__token)
        self.__code_string = self.__code_string[len(self.__token):]
        self.__code_string = self.__code_string.lstrip()
        self.__adjust_token()
        if log: self.log_token()

    def token(self) -> str:
        return self.__token

    def token_type(self) -> str:
        return self.__token_type

    def log_tokens(self) -> None:
        output_path = utils.init_output_file(self.__input_path, tokenizer_log=True)
        utils.write_in_file("<tokens>\n", output_path)
        while self.has_more_tokens():
            self.advance()
            utils.write_in_file("<" + self.__token_type + "> ", output_path)
            utils.write_in_file(self.__token, output_path)
            utils.write_in_file(" </" + self.__token_type + ">\n", output_path)
        utils.write_in_file("</tokens>", output_path)

    def log_token(self):
        utils.write_in_file("<" + self.__token_type + "> ", self.__output_path)
        utils.write_in_file(self.__token, self.__output_path)
        utils.write_in_file(" </" + self.__token_type + ">\n", self.__output_path)

    def __process_file(self) -> str:
        with open(self.__input_path, "r") as file:
            lines = file.readlines()

        file_content = []
        for line in lines:
            line = line.split("//", 1)[0]
            line = line.strip()
            if line:
                file_content.append(line)

        content_string = "".join(file_content)
        # remove /* type of comments
        while True:
            comment_start_index = content_string.find("/*")
            if comment_start_index == -1:
                break
            comment_end_index = content_string.find("*/", comment_start_index + 2)
            comment = content_string[comment_start_index: comment_end_index + 2]
            content_string = content_string.replace(comment, "")

        return content_string

    def __adjust_token(self):
        self.__token = self.__token.strip()
        if self.__token_type == "stringConstant":
            self.__token = self.__token.strip("\"")
        elif self.__token_type == "symbol":
            self.__token = self.__token.replace("&", "&amp;")
            self.__token = self.__token.replace("<", "&lt;")
            self.__token = self.__token.replace(">", "&gt;")
            self.__token = self.__token.replace("\"", "&quot;")
