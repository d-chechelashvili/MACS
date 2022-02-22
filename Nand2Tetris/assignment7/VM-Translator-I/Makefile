.PHONY: help
.DEFAULT_GOAL := help

help:
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

install: ## Install requirements
	pip install -r requirements.txt

fmt format: ## Run code formatters
	isort app vm_translator.py
	black app vm_translator.py

lint: ## Run code linters
	isort --check app vm_translator.py
	black --check app vm_translator.py
	flake8 app vm_translator.py
	mypy app vm_translator.py
