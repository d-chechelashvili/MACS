.PHONY: help
.DEFAULT_GOAL := help

help:
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

install: ## Install requirements
	pip install -r requirements.txt

fmt format: ## Run code formatters
	isort app tests spore.py
	black app tests spore.py

lint: ## Run code linters
	isort --check app tests spore.py
	black --check app tests spore.py
	flake8 app tests spore.py
	mypy app tests spore.py
