.PHONY: help
.DEFAULT_GOAL := help

help:
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

install: ## Install requirements
	pip install -r requirements.txt

fmt format: ## Run code formatters
	isort app tests point_of_sales.py
	black app tests point_of_sales.py

lint: ## Run code linters
	isort --check app tests point_of_sales.py
	black --check app tests point_of_sales.py
	flake8 app tests point_of_sales.py
	mypy app tests point_of_sales.py
