from app import utils


class VMWriter:
    def __init__(self, input_path):
        self.__output_path = utils.init_output_file(input_path)

    def write(self, text):
        utils.write_in_file(text, self.__output_path)

    def write_push(self, segment, index):
        self.write("push " + segment + " " + str(index) + "\n")

    def write_pop(self, segment, index):
        self.write("pop " + segment + " " + str(index) + "\n")

    def write_arithmetic(self, command):
        self.write(command + "\n")

    def write_label(self, label):
        self.write("label " + label + "\n")

    def write_go_to(self, label):
        self.write("goto " + label + "\n")

    def write_if(self, label):
        self.write("if-goto " + label + "\n")

    def write_call(self, name, args_num):
        self.write("call " + name + " " + str(args_num) + "\n")

    def write_function(self, name, local_var_num):
        self.write("function " + name + " " + str(local_var_num) + "\n")

    def write_return(self):
        self.write("return\n")