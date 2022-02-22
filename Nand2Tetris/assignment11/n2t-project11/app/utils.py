import os


def init_output_file(input_path: str) -> str:
    # output_path = os.path.join(os.getcwd(), "output")
    # if not os.path.exists(output_path):
    #     os.makedirs(output_path)
    output_path, file_name = os.path.split(input_path)
    output_name = file_name.split(".", 1)[0] + ".vm"
    output_path = os.path.join(output_path, output_name)
    with open(output_path, "w") as output_file:
        output_file.write("")
    return output_path


def write_in_file(content: str, output_path: str) -> None:
    with open(output_path, "a") as output_file:
        output_file.write(content)