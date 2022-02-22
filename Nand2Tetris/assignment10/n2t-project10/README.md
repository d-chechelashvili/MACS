## Conventions

Keep python conventions in mind:
  - variable_name
  - ClassName
  - CONSTANT_NAME
  - file_name.py

## Linters/Formatters

To avoid endless arguments related to coding style and formatting,
use "black" auto formatter and "isort" to order your imports.

Use mypy to check your static types. This not only helps you to catch errors but also
makes code better readable and understandable for your peers.

Use flake8 to catch semantic errors and common style issues.

Configuration for all the tools mentioned above is provided with this template project.

You can use `make` to run each of these tools or
see how to run them manually inside the `Makefile`.

## Requirements

Use following command to install needed requirements `pip install -r requirements.txt`
you can add your own development or production requirements as well.

## Usage

Use following command to see usage instructions `python jack_analyzer.py --help`
